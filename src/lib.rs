use wasm_bindgen::prelude::*;
use serde::{Serialize, Deserialize};
use serde_wasm_bindgen::to_value;
use wasm_bindgen::JsValue;

#[derive(Serialize, Deserialize)]
pub struct AlignmentResult {
    pub aligned_seq1: String,
    pub aligned_seq2: String,
    pub aligned_length: usize,
    pub aligned_identity: f64,
}

#[wasm_bindgen]
pub fn smith_waterman(seq1: &str, seq2: &str) -> JsValue {
    let result = smith_waterman_internal(seq1, seq2);
    to_value(&result).unwrap()
}

fn smith_waterman_internal(seq1: &str, seq2: &str) -> AlignmentResult {
    let seq1 = seq1.as_bytes();
    let seq2 = seq2.as_bytes();
    let len1 = seq1.len();
    let len2 = seq2.len();

    let match_score = 2;
    let mismatch_penalty = -1;
    let gap_penalty = -1;

    // Initialize scoring matrix
    let mut score_matrix = vec![vec![0; len2 + 1]; len1 + 1];

    let mut max_score = 0;
    let mut max_pos = (0, 0);

    // Fill the scoring matrix
    for i in 1..=len1 {
        for j in 1..=len2 {
            let match_mismatch = if seq1[i - 1] == seq2[j - 1] {
                match_score
            } else {
                mismatch_penalty
            };
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

    // Traceback
    let mut i = len1;
    let mut j = len2;
    let mut aligned_seq1 = Vec::new();
    let mut aligned_seq2 = Vec::new();
    let mut aligned_length = 0;
    let mut aligned_identity = 0.0;

    while i > max_pos.0  {
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

        if current_score == diagonal_score + if seq1[i - 1] == seq2[j - 1] { match_score } else { mismatch_penalty } {
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

