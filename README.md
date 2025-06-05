# unique_kmers_evolution

A Rust-based tool for visualizing the number of solid canonical kmers while streaming reads.

---

## Features

- Supports **FASTA** and **FASTQ**
- Supports **gzip-compressed** files (`.gz`)
- Real-time WebSocket output for monitoring
- **Growth**: The number of new solid k-mers between read intervals.
- **Acceleration**: The second derivative of k-mer discovery, indicating whether the rate of discovery is increasing, decreasing, or stabilizing
- Early termination based on configurable acceleration threshold


---

## Installation 

```bash

git clone https://github.com/yourusername/unique_kmers_evolution.git
cd unique_kmers_evolution
cargo build --release
cargo install --path .  
```

## Usage

```bash
unique_kmers_evolution --k 25 --input input.fq.gz
```


**Visualize** the evolution of the results opening file `plot.html` in a browser (reload the page once the program runs)
