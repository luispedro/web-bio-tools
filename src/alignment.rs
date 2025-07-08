use serde::{Serialize, Deserialize};

const BLOSUM62_MATRIX: [[i32; 20]; 20] = [
    [4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
    [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
    [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
    [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
    [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
    [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
    [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
    [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],
    [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
    [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],
    [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4],
];

fn aa_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'R' => Some(1),
        b'N' => Some(2),
        b'D' => Some(3),
        b'C' => Some(4),
        b'Q' => Some(5),
        b'E' => Some(6),
        b'G' => Some(7),
        b'H' => Some(8),
        b'I' => Some(9),
        b'L' => Some(10),
        b'K' => Some(11),
        b'M' => Some(12),
        b'F' => Some(13),
        b'P' => Some(14),
        b'S' => Some(15),
        b'T' => Some(16),
        b'W' => Some(17),
        b'Y' => Some(18),
        b'V' => Some(19),
        _ => None,
    }
}

fn blosum62_score(a: u8, b: u8) -> i32 {
    match (aa_index(a), aa_index(b)) {
        (Some(i), Some(j)) => BLOSUM62_MATRIX[i][j],
        _ => -4,
    }
}

#[derive(Serialize, Deserialize)]
pub struct AlignmentResult {
    pub aligned_seq1: String,
    pub aligned_seq2: String,
    pub aligned_length: usize,
    pub aligned_identity: f64,
}

pub fn smith_waterman_with_matrix<F>(
    seq1: &str,
    seq2: &str,
    gap_penalty: i32,
    score_fn: F,
) -> AlignmentResult
where
    F: Fn(u8, u8) -> i32,
{
    let seq1 = seq1.as_bytes();
    let seq2 = seq2.as_bytes();
    let len1 = seq1.len();
    let len2 = seq2.len();

    let mut score_matrix = vec![vec![0; len2 + 1]; len1 + 1];

    let mut max_score = 0;
    let mut max_pos = (0, 0);

    for i in 1..=len1 {
        for j in 1..=len2 {
            let match_mismatch = score_fn(seq1[i - 1], seq2[j - 1]);
            let diagonal_score = score_matrix[i - 1][j - 1] + match_mismatch;
            let up_score = score_matrix[i - 1][j] + gap_penalty;
            let left_score = score_matrix[i][j - 1] + gap_penalty;
            let cell_score = diagonal_score.max(up_score).max(left_score).max(0);
            score_matrix[i][j] = cell_score;

            if cell_score > max_score {
                max_score = cell_score;
                max_pos = (i, j);
            }
        }
    }

    let mut i = len1;
    let mut j = len2;
    let mut aligned_seq1 = Vec::new();
    let mut aligned_seq2 = Vec::new();
    let mut aligned_length = 0;
    let mut aligned_identity = 0.0;

    while i > max_pos.0 {
        aligned_seq1.push(seq1[i - 1]);
        aligned_seq2.push(b'-');
        i -= 1;
    }

    while j > max_pos.1 {
        aligned_seq1.push(b'-');
        aligned_seq2.push(seq2[j - 1]);
        j -= 1;
    }

    while i > 0 && j > 0 && score_matrix[i][j] > 0 {
        let current_score = score_matrix[i][j];
        let diagonal_score = score_matrix[i - 1][j - 1];
        let up_score = score_matrix[i - 1][j];
        let left_score = score_matrix[i][j - 1];

        if current_score == diagonal_score + score_fn(seq1[i - 1], seq2[j - 1]) {
            aligned_seq1.push(seq1[i - 1]);
            aligned_seq2.push(seq2[j - 1]);
            i -= 1;
            j -= 1;
            aligned_length += 1;
            if seq1[i] == seq2[j] {
                aligned_identity += 1.0;
            }
        } else if current_score == up_score + gap_penalty {
            aligned_seq1.push(seq1[i - 1]);
            aligned_seq2.push(b'-');
            i -= 1;
        } else if current_score == left_score + gap_penalty {
            aligned_seq1.push(b'-');
            aligned_seq2.push(seq2[j - 1]);
            j -= 1;
        } else {
            break;
        }
    }
    while i > 0 {
        aligned_seq1.push(seq1[i - 1]);
        aligned_seq2.push(b'-');
        i -= 1;
    }
    while j > 0 {
        aligned_seq1.push(b'-');
        aligned_seq2.push(seq2[j - 1]);
        j -= 1;
    }

    aligned_seq1.reverse();
    aligned_seq2.reverse();

    AlignmentResult {
        aligned_seq1: String::from_utf8(aligned_seq1).unwrap(),
        aligned_seq2: String::from_utf8(aligned_seq2).unwrap(),
        aligned_length,
        aligned_identity: if aligned_length > 0 {
            aligned_identity / aligned_length as f64
        } else {
            0.0
        },
    }
}

pub fn smith_waterman_blosum62_internal(
    seq1: &str,
    seq2: &str,
    gap_penalty: i32,
) -> AlignmentResult {
    smith_waterman_with_matrix(seq1, seq2, gap_penalty, blosum62_score)
}

pub fn smith_waterman_internal(
    seq1: &str,
    seq2: &str,
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> AlignmentResult {
    smith_waterman_with_matrix(seq1, seq2, gap_penalty, |a, b| {
        if a == b { match_score } else { mismatch_penalty }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn align_identical_sequences() {
        let r = smith_waterman_internal("GATTACA", "GATTACA", 2, -1, -1);
        assert_eq!(r.aligned_seq1, "GATTACA");
        assert_eq!(r.aligned_seq2, "GATTACA");
        assert_eq!(r.aligned_length, 7);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn align_acacacta_agcacaca() {
        let r = smith_waterman_internal("ACACACTA", "AGCACACA", 2, -1, -1);
        assert_eq!(r.aligned_seq1, "A-CACACTA");
        assert_eq!(r.aligned_seq2, "AGCACAC-A");
        assert_eq!(r.aligned_length, 7);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn align_gattaca_gcatgcu() {
        let r = smith_waterman_internal("GATTACA", "GCATGCU", 2, -1, -1);
        assert_eq!(r.aligned_seq1, "G-AT---TACA");
        assert_eq!(r.aligned_seq2, "GCATGCU----");
        assert_eq!(r.aligned_length, 3);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn blosum62_identical_sequences() {
        let r = smith_waterman_blosum62_internal("GATTACA", "GATTACA", -1);
        assert_eq!(r.aligned_seq1, "GATTACA");
        assert_eq!(r.aligned_seq2, "GATTACA");
        assert_eq!(r.aligned_length, 7);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
    }
}
