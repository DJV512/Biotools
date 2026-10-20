use crate::constants::BLOSUMX;
use crate::constants::aa_index;

#[derive(Clone, Copy, Debug, PartialEq)]
enum Direction {
    Diag,
    Vert,
    Hori,
    Start,
    Impossible,
}


pub fn pairwise_protein(seq1: &str, seq2: &str) -> (String, String) {
    // println!("Aligning {seq1} against {seq2}");

    // Set Gap penalty constants, and a highly negative number for the impossible boundaries
    const GAP_OPEN:i32 = -10;
    const GAP_EXTEND:i32 = -1;
    const NEG_INF:i32 = i32::MIN / 2;

    // Add a leading space to the sequences for scoring
    let mut s1 = " ".to_string();
    let mut s2 = " ".to_string();
    s1.push_str(seq1);
    s2.push_str(seq2);

    // Determine sequence lengths and the necessary size of the arrays
    let rows = s1.len();
    let cols = s2.len();
    let size = rows * cols;

    // Initialize scoring and traceback matrices
    let mut score_diag:Vec<i32> = vec![0; size];
    let mut score_vert:Vec<i32> = vec![0; size];
    let mut score_hori:Vec<i32> = vec![0; size];
    let mut trace_diag: Vec<Direction> = vec![Direction::Impossible; size];
    let mut trace_vert: Vec<Direction> = vec![Direction::Impossible; size];
    let mut trace_hori: Vec<Direction> = vec![Direction::Impossible; size];

    
    // Set initial values for the array origins
    score_vert[0] = NEG_INF;
    score_hori[0] = NEG_INF;
    trace_diag[0] = Direction::Start;
    trace_vert[0] = Direction::Impossible;
    trace_hori[0] = Direction::Impossible;


    // Set initial value for the array top rows
    for j in 1..cols {
        score_diag[j] = NEG_INF;
        score_vert[j] = NEG_INF;
        score_hori[j] = GAP_OPEN + (j as i32 - 1) * GAP_EXTEND;

        trace_diag[j] = Direction::Impossible;
        trace_vert[j] = Direction::Impossible;
        
        if j == 1 {
            trace_hori[j] = Direction::Diag;
        } else {
            trace_hori[j] = Direction::Hori;
        }
        
    }

    // Set initial values for the array first columns
    for i in 1..rows {
        score_diag[i*cols] = NEG_INF;
        score_vert[i*cols] = GAP_OPEN + (i as i32 - 1) * GAP_EXTEND;
        score_hori[i*cols] = NEG_INF;

        trace_diag[i*cols] = Direction::Impossible;
        trace_hori[i*cols] = Direction::Impossible;

        if i == 1 {
            trace_vert[i*cols] = Direction::Diag;
        } else {
            trace_vert[i*cols] = Direction::Vert;
        }
    }

    // Populate the entire scoring and traceback matrix based on maximum BLOSUM62 based scoring
    // BLOSUMX is a modified BLOSUM62 whcih uses an "X" amino acid as an alignment anchor
    for i in 1..rows {
        for j in 1..cols {
            let residue_blosum_score = BLOSUMX[aa_index(&s1[i..i+1])][aa_index(&s2[j..j+1])];

            let (best_score, direction) = [
                (score_vert[(i-1)*cols + (j-1)] + residue_blosum_score, Direction::Vert),
                (score_hori[(i-1)*cols + (j-1)] + residue_blosum_score, Direction::Hori),
                (score_diag[(i-1)*cols + (j-1)] + residue_blosum_score, Direction::Diag), 
                ].into_iter().max_by_key(|(score, _)| *score).unwrap();

            score_diag[i*cols + j] = best_score;
            trace_diag[i*cols + j] = direction;

            let (best_score, direction) = [
                (score_vert[(i-1)*cols + j] + GAP_EXTEND, Direction::Vert),
                (score_diag[(i-1)*cols + j] + GAP_OPEN, Direction::Diag),
            ].into_iter().max_by_key(|(score, _)| *score).unwrap();

            score_vert[i*cols + j] = best_score;
            trace_vert[i*cols + j] = direction;

            let (best_score, direction) = [
                (score_hori[i*cols + (j-1)] + GAP_EXTEND, Direction::Hori),
                (score_diag[i*cols + (j-1)] + GAP_OPEN, Direction::Diag),
            ].into_iter().max_by_key(|(score, _)| *score).unwrap();

            score_hori[i*cols + j] = best_score;
            trace_hori[i*cols + j] = direction;
        }
    }

    // Determine the best score at the alignment end from the three scoring matrices
    let last = size - 1;
    let (_, traceback_matrix_start) = [
        (score_vert[last], Direction::Vert),
        (score_hori[last], Direction::Hori),
        (score_diag[last], Direction::Diag),
    ].into_iter().max_by_key(|(score, _)| *score).unwrap();

    traceback(&s1, &s2, &trace_diag, &trace_hori, &trace_vert, traceback_matrix_start, rows, cols)
}


fn traceback(
    s1:&str, 
    s2:&str, 
    trace_diag:&Vec<Direction>, 
    trace_hori:&Vec<Direction>, 
    trace_vert:&Vec<Direction>, 
    traceback_matrix_start:Direction, 
    rows:usize, 
    cols:usize
) -> (String, String) {

    let mut align1 = "".to_string();
    let mut align2 = "".to_string();
    let mut predecessor = traceback_matrix_start;
    let mut i = rows-1;
    let mut j = cols-1;

    while i > 0 || j > 0 {

        match predecessor {
            Direction::Diag => {
                predecessor = trace_diag[i*cols + j];
                align1.push_str(&s1[i..i+1]);
                align2.push_str(&s2[j..j+1]);
                i -= 1;
                j -= 1;
            }
            Direction::Hori => {
                predecessor = trace_hori[i*cols + j];
                align1.push('-');
                align2.push_str(&s2[j..j+1]);
                j -= 1;
            }
            Direction::Vert => {
                predecessor = trace_vert[i*cols + j];
                align1.push_str(&s1[i..i+1]);
                align2.push('-');
                i -= 1;
            }
            Direction::Start => {
                break
            }
            Direction::Impossible => {
                panic!("Traceback error: impossible state at ({i}, {j})")
            }
        }
    }

    align1 = align1.chars().rev().collect();
    align2 = align2.chars().rev().collect();

    assert_eq!(align1.len(), align2.len());

    (align1, align2)
}