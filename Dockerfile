# docker build -t unique_kmers_evolution .
# docker run --rm unique_kmers_evolution

FROM rust:1.87

WORKDIR /usr/src/unique_kmers_evolution
COPY . .

RUN cargo build --release
RUN cargo install --path .

ENTRYPOINT ["unique_kmers_evolution"]
