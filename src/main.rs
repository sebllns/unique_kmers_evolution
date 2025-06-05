use clap::Parser;
use futures::{SinkExt, StreamExt};
use std::path::PathBuf;
use std::sync::Arc;
use tokio::sync::{mpsc, Mutex};
use warp::ws::{Message, WebSocket};
use warp::Filter;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use bio::io::{fasta, fastq};


/// Fast hash map
use rustc_hash::FxHashMap;

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Length of k-mers
    #[arg(short, long)]
    k: usize,

    /// Input FASTA file
    #[arg(short, long)]
    input: PathBuf,
}

enum RecordReader<R: Read> {
    Fasta(fasta::Records<BufReader<R>>),
    Fastq(fastq::Records<BufReader<R>>),
}

impl<R: Read> RecordReader<R> {
    fn next_record(&mut self) -> Option<Result<Vec<u8>, Box<dyn std::error::Error>>> {
        match self {
            RecordReader::Fasta(reader) => reader.next().map(|r| {
                r.map(|rec| rec.seq().to_vec())
                    .map_err(|e| e.into())
            }),
            RecordReader::Fastq(reader) => reader.next().map(|r| {
                r.map(|rec| rec.seq().to_vec())
                    .map_err(|e| e.into())
            }),
        }
    }
}

/// Fast reverse complement for &[u8]
fn reverse_complement(kmer: &[u8]) -> Vec<u8> {
    kmer.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => c,
        })
        .collect()
}

/// Return the canonical k-mer (lexicographically smallest between kmer and its reverse complement)
fn canonical_kmer(kmer: &[u8]) -> Vec<u8> {
    let rc = reverse_complement(kmer);
    if kmer <= rc.as_slice() {
        kmer.to_vec()
    } else {
        rc
    }
}

/// WebSocket handling
async fn handle_connection(ws: WebSocket, rx: Arc<Mutex<mpsc::Receiver<(u32, u32)>>>) {
    let (mut ws_tx, _) = ws.split();
    let mut rx = rx.lock().await;
    while let Some((reads, kmers)) = rx.recv().await {
        let message = format!("{} {}", reads, kmers);
        if ws_tx.send(Message::text(message)).await.is_err() {
            break;
        }
    }
}


fn open_reader(path: &PathBuf) -> Result<RecordReader<impl Read>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader: Box<dyn Read> = if path.extension().map(|e| e == "gz").unwrap_or(false) {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut buffered = BufReader::new(reader);

    // Peek at the first byte
    let first_byte = {
        let buf = buffered.fill_buf()?;
        if buf.is_empty() {
            return Err("Input file is empty".into());
        }
        buf[0]
    };

    // Decide format by first byte
    if first_byte == b'>' {
        Ok(RecordReader::Fasta(fasta::Reader::new(buffered).records()))
    } else if first_byte == b'@' {
        Ok(RecordReader::Fastq(fastq::Reader::new(buffered).records()))
    } else {
        Err(format!("Unknown file format: expected '>' or '@', got '{}'", first_byte as char).into())
    }
}


#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let k = args.k;

    let mut unique_kmers: FxHashMap<Vec<u8>, bool> = FxHashMap::default();
    let mut unique_solid_kmers = 0;

    let (tx, rx) = mpsc::channel(100);
    let rx = Arc::new(Mutex::new(rx));

    // WebSocket route
    let ws_route = warp::path("ws")
        .and(warp::ws())
        .map(move |ws: warp::ws::Ws| {
            let rx = rx.clone();
            ws.on_upgrade(move |socket| handle_connection(socket, rx))
        });

    tokio::spawn(async move {
        warp::serve(ws_route)
            .run(([127, 0, 0, 1], 3030))
            .await;
    });

    let mut reader = open_reader(&args.input)?;
    let mut idx = 0;

    let mut prev_kmers = 0u32;
    let mut growth_history: Vec<i32> = Vec::new();
    let mut accel_history: Vec<i32> = Vec::new();

    while let Some(seq_result) = reader.next_record() {
        let sequence = seq_result?;

        if sequence.len() < k {
            idx += 1;
            continue;
        }

        for i in 0..=(sequence.len() - k) {
            let kmer = &sequence[i..i + k];
            let canonical = canonical_kmer(kmer);
            match unique_kmers.get_mut(&canonical) {
                Some(seen) => {
                    if !*seen {
                        *seen = true;
                        unique_solid_kmers += 1;
                    }
                }
                None => {
                    unique_kmers.insert(canonical, false);
                }
            }
        }

        if idx % 10000 == 0 {

            let reads = idx as u32;
            let kmers = unique_solid_kmers;
            let growth = kmers as i32 - prev_kmers as i32;

            growth_history.push(growth);
            if growth_history.len() > 10 {
                growth_history.remove(0);
            }

            // Compute acceleration only if we have at least 2 growth values
            if growth_history.len() >= 2 {
                let acceleration = growth_history[growth_history.len() - 1]
                    - growth_history[growth_history.len() - 2];
                accel_history.push(acceleration);
                if accel_history.len() > 10 {
                    accel_history.remove(0);
                }
            }

            // Compute averages
            let avg_growth = growth_history.iter().sum::<i32>() as f32 / growth_history.len() as f32;
            let avg_accel = if !accel_history.is_empty() {
                accel_history.iter().sum::<i32>() as f32 / accel_history.len() as f32
            } else {
                0.0
            };

            println!(
                "Processed {} reads, unique k-mers: {}, Δ_avg: {:.1}, Δ²_avg: {:.1}",
                reads, kmers, avg_growth, avg_accel
            );

            // WebSocket message can include acceleration too if desired
            tx.send((reads, kmers)).await?;

            // Auto-stop condition
            if reads > 50000 && avg_accel.abs() < 20.0 {
                println!(
                    "Stopping early: acceleration average {:.1} < 50 after {} reads.",
                    avg_accel, reads
                );
                break;
            }

            prev_kmers = kmers;
        }

        idx += 1;
    }

    Ok(())
}
