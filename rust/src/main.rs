use biotools_rust::pairwise_protein;

fn main() {
    let seq1 = "AVMYERF";
    let seq2 = "QAVMFDKF";
    pairwise_protein(seq1, seq2);
}