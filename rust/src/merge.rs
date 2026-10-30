use crate::constants::{BLOSUMX,aa_index,Direction,Profile};


pub fn merge(left: &Profile, right: &Profile, gap_open: f64, gap_extend: f64) -> Profile {

    const NEG_INF:f64 = f64::NEG_INFINITY;  

    // Add a leading space to the sequences for scoring
    let mut _left = left.clone();    
    let mut _right = right.clone();
    for seq in &mut _left.sequences {
        seq.insert(0, '-');
    }
    for seq in &mut _right.sequences {
        seq.insert(0, '-');
    }

    // Determine sequence lengths and the necessary size of the arrays
    let rows = _left.sequences[0].len();
    let cols = _right.sequences[0].len();
    let size = rows * cols;

    // Initialize scoring and traceback matrices
    let mut score_diag:Vec<f64> = vec![0.0; size];
    let mut score_vert:Vec<f64> = vec![0.0; size];
    let mut score_hori:Vec<f64> = vec![0.0; size];
    let mut trace_diag: Vec<Direction> = vec![Direction::Impossible; size];
    let mut trace_vert: Vec<Direction> = vec![Direction::Impossible; size];
    let mut trace_hori: Vec<Direction> = vec![Direction::Impossible; size];

    // Set initial values for the array origins
    score_vert[0] = NEG_INF;
    score_hori[0] = NEG_INF;
    trace_diag[0] = Direction::Start;
    trace_vert[0] = Direction::Impossible;
    trace_hori[0] = Direction::Impossible;

    // Top row scores
    for j in 1..cols {
        score_diag[j] = NEG_INF;
        score_vert[j] = NEG_INF;
        score_hori[j] = gap_open + (j - 1) as f64 * gap_extend;
        trace_diag[j] = Direction::Impossible;
        trace_vert[j] = Direction::Impossible;
        if j == 1 {
            trace_hori[j] = Direction::Diag;
        } else {
            trace_hori[j] = Direction::Hori;
        }
    }

    // First column scores
    for i in 1..rows {
        score_diag[i * cols] = NEG_INF;
        score_hori[i * cols] = NEG_INF;
        score_vert[i * cols] = gap_open + (i - 1) as f64 * gap_extend;
        trace_diag[i * cols] = Direction::Impossible;
        trace_hori[i * cols] = Direction::Impossible;
        if i == 1 {
            trace_vert[i * cols] = Direction::Diag;
        } else {
            trace_vert[i * cols] = Direction::Vert;
        }
    }

    // Interior scores
    // Iterate of the length of profile 1
    for i in 1..rows {
        let left_column: Vec<char> = _left.sequences
            .iter()
            .map(|seq| seq[i])
            .collect();

        for j in 1..cols {
            let right_column: Vec<char> = _right.sequences
            .iter()
            .map(|seq| seq[j])
            .collect();


            // Current pair is a residue column vs residue column
            let column_v_column_blosum_score = _msa_sigma(&left_column, &right_column);

            let (best_score, direction) = [
                (score_vert[(i-1) * cols + (j-1)] + column_v_column_blosum_score, Direction::Vert),
                (score_hori[(i-1) * cols + (j-1)] + column_v_column_blosum_score, Direction::Hori),
                (score_diag[(i-1) * cols + (j-1)] + column_v_column_blosum_score, Direction::Diag),
            ].into_iter().max_by(|a, b| a.0.total_cmp(&b.0)).unwrap();
            score_diag[i * cols + j] = best_score;
            trace_diag[i * cols + j] = direction;

            // Current pair is a gap column in profile 2
            let (best_score, direction) = [
                (score_vert[(i-1)*cols + j] + gap_extend, Direction::Vert),
                (score_diag[(i-1)*cols + j] + gap_open, Direction::Diag),
            ].into_iter().max_by(|a, b| a.0.total_cmp(&b.0)).unwrap();
            score_vert[i * cols + j] = best_score;
            trace_vert[i * cols + j] = direction;

            // Current pair is a gap column in profile 1
            let (best_score, direction) = [
                (score_hori[(i * cols) + (j-1)] + gap_extend, Direction::Hori),
                (score_diag[(i * cols) + (j-1)] + gap_open, Direction::Diag),
            ].into_iter().max_by(|a, b| a.0.total_cmp(&b.0)).unwrap();
            score_hori[i * cols + j] = best_score;
            trace_hori[i * cols + j] = direction

        }
    }

    let last = size - 1;
    let (_, traceback_matrix_start) = [
        (score_hori[last], Direction::Hori),
        (score_vert[last], Direction::Vert),
        (score_diag[last], Direction::Diag),
    ].into_iter().max_by(|a, b| a.0.total_cmp(&b.0)).unwrap();

    msa_traceback(&_left, &_right, &trace_diag, &trace_hori, &trace_vert, traceback_matrix_start, rows, cols)
}


fn _msa_sigma(
    left_residues: &[char],
    right_residues: &[char],
) -> f64 {
    let comparisons = left_residues.len() * right_residues.len();
    let mut score = 0;
    for &left_residue in left_residues {
        for &right_residue in right_residues {
            if left_residue == '-' && right_residue == '-' {
                continue;
            } else if left_residue == '-' || right_residue == '-' {
                score -= 8;
            } else {
                let mut left_buf = [0u8; 4];
                let mut right_buf = [0u8; 4];   

                let left_aa = left_residue.encode_utf8(&mut left_buf);
                let right_aa = right_residue.encode_utf8(&mut right_buf);

                score += BLOSUMX[aa_index(left_aa)][aa_index(right_aa)];
            }
        }
    }
    score as f64 / comparisons as f64
}



fn msa_traceback(
    _left:&Profile,
    _right:&Profile,
    trace_diag:&Vec<Direction>, 
    trace_hori:&Vec<Direction>, 
    trace_vert:&Vec<Direction>,
    traceback_matrix_start:Direction,
    rows:usize, 
    cols:usize, 
) -> Profile {

    let total_sequences = _left.sequences.len() + _right.sequences.len();
    let total_left = _left.sequences.len();

    let names = _left.names
    .iter()
    .chain(_right.names.iter())
    .cloned()
    .collect();

    let mut merged_profile = Profile {
        names,
        sequences: vec![Vec::new(); total_sequences],
    };

    let mut predecessor = traceback_matrix_start;
    let mut i = rows-1;
    let mut j = cols-1;

    while i > 0 || j > 0 {

        match predecessor {
            Direction::Diag => {
                predecessor = trace_diag[i * cols + j];
                for (z, sequence) in _left.sequences.iter().enumerate() {
                    merged_profile.sequences[z].push(sequence[i]);
                }
                for (z, sequence) in _right.sequences.iter().enumerate() {
                    merged_profile.sequences[total_left + z].push(sequence[j]);
                }
                i -= 1;
                j -= 1;
            }
            Direction::Hori => {
                predecessor = trace_hori[i * cols + j];
                for (z, _) in _left.sequences.iter().enumerate() {
                    merged_profile.sequences[z].push('-');
                }
                for (z, sequence) in _right.sequences.iter().enumerate() {
                    merged_profile.sequences[total_left + z].push(sequence[j]);
                }
                j -= 1;
            }
            Direction::Vert => {
                predecessor = trace_vert[i * cols + j];
                for (z, sequence) in _left.sequences.iter().enumerate() {
                    merged_profile.sequences[z].push(sequence[i]);
                }
                for (z, _) in _right.sequences.iter().enumerate() {
                    merged_profile.sequences[total_left + z].push('-');
                }
                i -= 1;
            }
            Direction::Start => {
                break;
            }
            Direction::Impossible => {
                panic!("Traceback error: impossible state at ({i}, {j})")
            }    
        }
    }

    for sequence in &mut merged_profile.sequences {
        sequence.reverse();
    }

    merged_profile
}