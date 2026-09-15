use crate::alignment::{self, AlignmentResult};
use serde::{Deserialize, Serialize};

/// Scoring configuration shared by every pair of an all-by-all run.
#[derive(Clone, Copy)]
pub struct ScoringParams {
    pub global: bool,
    pub blosum62: bool,
    pub match_score: f64,
    pub mismatch_penalty: f64,
    pub gap_open: f64,
    pub gap_extend: f64,
}

/// Summary of one pairwise alignment, without the aligned sequences themselves.
#[derive(Serialize, Deserialize, Clone, PartialEq, Debug)]
pub struct PairSummary {
    pub identity: f64,
    pub score: f64,
    pub aligned_length: usize,
}

pub fn align_pair(seq1: &str, seq2: &str, params: &ScoringParams) -> AlignmentResult {
    match (params.global, params.blosum62) {
        (true, true) => alignment::needleman_wunsch_blosum62_internal(
            seq1,
            seq2,
            params.gap_open,
            params.gap_extend,
        ),
        (true, false) => alignment::needleman_wunsch_internal(
            seq1,
            seq2,
            params.match_score,
            params.mismatch_penalty,
            params.gap_open,
            params.gap_extend,
        ),
        (false, true) => alignment::smith_waterman_blosum62_internal(
            seq1,
            seq2,
            params.gap_open,
            params.gap_extend,
        ),
        (false, false) => alignment::smith_waterman_internal(
            seq1,
            seq2,
            params.match_score,
            params.mismatch_penalty,
            params.gap_open,
            params.gap_extend,
        ),
    }
}

/// Aligns `seqs[row]` against itself and every sequence after it.
///
/// Entry `k` of the result is the alignment of `seqs[row]` with `seqs[row + k]`,
/// so the first entry is always the self-alignment. Only this upper triangle is
/// computed: the caller mirrors it into the lower triangle. Alignment scores are
/// symmetric, and where a tie in the traceback would let the transposed
/// alignment report a different identity, either answer is equally good, so a
/// mirrored matrix is preferable to an asymmetric one.
pub fn all_by_all_row(seqs: &[String], row: usize, params: &ScoringParams) -> Vec<PairSummary> {
    seqs[row..]
        .iter()
        .map(|other| {
            let result = align_pair(&seqs[row], other, params);
            PairSummary {
                identity: result.aligned_identity,
                score: result.score,
                aligned_length: result.aligned_length,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blosum_params(global: bool) -> ScoringParams {
        ScoringParams {
            global,
            blosum62: true,
            match_score: 2.0,
            mismatch_penalty: -1.0,
            gap_open: -10.0,
            gap_extend: -0.5,
        }
    }

    fn seqs() -> Vec<String> {
        vec![
            "MKTAYIAKQRQISFVKSHFSRQ".to_string(),
            "MKTAYIAKQRQISFVKSHFSRQLEER".to_string(),
            "WWWWPPPPCCCC".to_string(),
        ]
    }

    #[test]
    fn row_covers_upper_triangle_including_diagonal() {
        let seqs = seqs();
        let params = blosum_params(false);
        assert_eq!(all_by_all_row(&seqs, 0, &params).len(), 3);
        assert_eq!(all_by_all_row(&seqs, 1, &params).len(), 2);
        assert_eq!(all_by_all_row(&seqs, 2, &params).len(), 1);
    }

    #[test]
    fn self_alignment_is_fully_identical() {
        let seqs = seqs();
        for global in [false, true] {
            let params = blosum_params(global);
            for row in 0..seqs.len() {
                let diagonal = &all_by_all_row(&seqs, row, &params)[0];
                assert!((diagonal.identity - 1.0).abs() < 1e-9);
                assert_eq!(diagonal.aligned_length, seqs[row].len());
            }
        }
    }

    #[test]
    fn scores_are_symmetric() {
        let seqs = seqs();
        for global in [false, true] {
            let params = blosum_params(global);
            for i in 0..seqs.len() {
                let row = all_by_all_row(&seqs, i, &params);
                for (k, entry) in row.iter().enumerate() {
                    let transposed = align_pair(&seqs[i + k], &seqs[i], &params);
                    assert!((entry.score - transposed.score).abs() < 1e-6);
                }
            }
        }
    }

    #[test]
    fn similar_sequences_score_above_unrelated_ones() {
        let seqs = seqs();
        let row = all_by_all_row(&seqs, 0, &blosum_params(false));
        assert!(row[1].score > row[2].score);
        assert!(row[1].identity > row[2].identity);
    }

    #[test]
    fn row_matches_the_pairwise_functions() {
        let seqs = seqs();
        let params = blosum_params(false);
        let row = all_by_all_row(&seqs, 0, &params);
        let direct = alignment::smith_waterman_blosum62_internal(
            &seqs[0],
            &seqs[1],
            params.gap_open,
            params.gap_extend,
        );
        assert!((row[1].identity - direct.aligned_identity).abs() < 1e-9);
        assert!((row[1].score - direct.score).abs() < 1e-9);
        assert_eq!(row[1].aligned_length, direct.aligned_length);
    }

    #[test]
    fn uniform_scoring_is_used_when_blosum62_is_off() {
        let dna = vec!["ACGTACGT".to_string(), "ACGTACGT".to_string()];
        let params = ScoringParams {
            global: false,
            blosum62: false,
            match_score: 2.0,
            mismatch_penalty: -1.0,
            gap_open: -10.0,
            gap_extend: -0.5,
        };
        let row = all_by_all_row(&dna, 0, &params);
        assert!((row[1].score - 16.0).abs() < 1e-9);
    }
}
