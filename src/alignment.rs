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

fn blosum62_score(a: u8, b: u8) -> f64 {
    match (aa_index(a), aa_index(b)) {
        (Some(i), Some(j)) => BLOSUM62_MATRIX[i][j] as f64,
        _ => -4.0,
    }
}

#[derive(Serialize, Deserialize, Clone)]
pub struct AlignmentResult {
    pub aligned_seq1: String,
    pub aligned_seq2: String,
    pub aligned_length: usize,
    pub aligned_identity: f64,
    pub score: f64,
    pub alignment_markup: String,
}

pub fn smith_waterman_with_matrix<F>(
    seq1: &str,
    seq2: &str,
    gap_open: f64,
    gap_extend: f64,
    score_fn: F,
) -> AlignmentResult
where
    F: Fn(u8, u8) -> f64,
{
    let seq1 = seq1.as_bytes();
    let seq2 = seq2.as_bytes();
    let len1 = seq1.len();
    let len2 = seq2.len();

    let mut score_matrix = vec![vec![0.0; len2 + 1]; len1 + 1];
    let mut ins_matrix = vec![vec![f64::NEG_INFINITY; len2 + 1]; len1 + 1];
    let mut del_matrix = vec![vec![f64::NEG_INFINITY; len2 + 1]; len1 + 1];

    let mut max_score = 0.0;
    let mut max_pos = (0, 0);

    for i in 1..=len1 {
        for j in 1..=len2 {
            let match_mismatch = score_fn(seq1[i - 1], seq2[j - 1]);
            let diagonal_score = score_matrix[i - 1][j - 1] + match_mismatch;

            ins_matrix[i][j] =
                (score_matrix[i - 1][j] + gap_open).max(ins_matrix[i - 1][j] + gap_extend);
            del_matrix[i][j] =
                (score_matrix[i][j - 1] + gap_open).max(del_matrix[i][j - 1] + gap_extend);

            let cell_score = diagonal_score
                .max(ins_matrix[i][j])
                .max(del_matrix[i][j])
                .max(0.0);
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

    while i > 0 && j > 0 && score_matrix[i][j] > 0.0 {
        let current_score = score_matrix[i][j];
        let diagonal_score = score_matrix[i - 1][j - 1];
        let up_score = ins_matrix[i][j];
        let left_score = del_matrix[i][j];

        let diag = diagonal_score + score_fn(seq1[i - 1], seq2[j - 1]);
        if (current_score - diag).abs() < 1e-6 {
            aligned_seq1.push(seq1[i - 1]);
            aligned_seq2.push(seq2[j - 1]);
            i -= 1;
            j -= 1;
            aligned_length += 1;
            if seq1[i].to_ascii_uppercase() == seq2[j].to_ascii_uppercase() {
                aligned_identity += 1.0;
            }
        } else if (current_score - up_score).abs() < 1e-6 {
            aligned_seq1.push(seq1[i - 1]);
            aligned_seq2.push(b'-');
            i -= 1;
        } else if (current_score - left_score).abs() < 1e-6 {
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

    let mut markup = String::with_capacity(aligned_seq1.len());
    for (&a, &b) in aligned_seq1.iter().zip(aligned_seq2.iter()) {
        let ch = if a.to_ascii_uppercase() == b.to_ascii_uppercase() && a != b'-' && b != b'-' {
            '|'
        } else if a == b'-' || b == b'-' {
            ' '
        } else if score_fn(a, b) > 0.0 {
            ':'
        } else {
            '.'
        };
        markup.push(ch);
    }

    AlignmentResult {
        aligned_seq1: String::from_utf8(aligned_seq1).unwrap(),
        aligned_seq2: String::from_utf8(aligned_seq2).unwrap(),
        aligned_length,
        aligned_identity: if aligned_length > 0 {
            aligned_identity / aligned_length as f64
        } else {
            0.0
        },
        score: max_score,
        alignment_markup: markup,
    }
}

pub fn smith_waterman_blosum62_internal(
    seq1: &str,
    seq2: &str,
    gap_open: f64,
    gap_extend: f64,
) -> AlignmentResult {
    smith_waterman_with_matrix(seq1, seq2, gap_open, gap_extend, |a, b| {
        blosum62_score(a, b) as f64
    })
}

pub fn smith_waterman_internal(
    seq1: &str,
    seq2: &str,
    match_score: f64,
    mismatch_penalty: f64,
    gap_open: f64,
    gap_extend: f64,
) -> AlignmentResult {
    smith_waterman_with_matrix(seq1, seq2, gap_open, gap_extend, |a, b| {
        if a.to_ascii_uppercase() == b.to_ascii_uppercase() {
            match_score
        } else {
            mismatch_penalty
        }
    })
}

pub fn needleman_wunsch_with_matrix<F>(
    seq1: &str,
    seq2: &str,
    gap_open: f64,
    gap_extend: f64,
    score_fn: F,
) -> AlignmentResult
where
    F: Fn(u8, u8) -> f64,
{
    let seq1 = seq1.as_bytes();
    let seq2 = seq2.as_bytes();
    let len1 = seq1.len();
    let len2 = seq2.len();

    let mut score_matrix = vec![vec![0.0; len2 + 1]; len1 + 1];
    let mut ins_matrix = vec![vec![f64::NEG_INFINITY; len2 + 1]; len1 + 1];
    let mut del_matrix = vec![vec![f64::NEG_INFINITY; len2 + 1]; len1 + 1];

    for i in 1..=len1 {
        ins_matrix[i][0] = if i == 1 {
            score_matrix[i - 1][0] + gap_open
        } else {
            ins_matrix[i - 1][0] + gap_extend
        };
        score_matrix[i][0] = ins_matrix[i][0];
    }
    for j in 1..=len2 {
        del_matrix[0][j] = if j == 1 {
            score_matrix[0][j - 1] + gap_open
        } else {
            del_matrix[0][j - 1] + gap_extend
        };
        score_matrix[0][j] = del_matrix[0][j];
    }

    for i in 1..=len1 {
        for j in 1..=len2 {
            let match_mismatch = score_fn(seq1[i - 1], seq2[j - 1]);
            let diagonal_score = score_matrix[i - 1][j - 1] + match_mismatch;

            ins_matrix[i][j] =
                (score_matrix[i - 1][j] + gap_open).max(ins_matrix[i - 1][j] + gap_extend);
            del_matrix[i][j] =
                (score_matrix[i][j - 1] + gap_open).max(del_matrix[i][j - 1] + gap_extend);

            score_matrix[i][j] = diagonal_score.max(ins_matrix[i][j]).max(del_matrix[i][j]);
        }
    }

    let mut i = len1;
    let mut j = len2;
    let mut aligned_seq1 = Vec::new();
    let mut aligned_seq2 = Vec::new();
    let mut aligned_length = 0;
    let mut aligned_identity = 0.0;

    while i > 0 || j > 0 {
        if i > 0
            && j > 0
            && (score_matrix[i][j]
                - (score_matrix[i - 1][j - 1] + score_fn(seq1[i - 1], seq2[j - 1])))
                .abs()
                < 1e-6
        {
            aligned_seq1.push(seq1[i - 1]);
            aligned_seq2.push(seq2[j - 1]);
            if seq1[i - 1].to_ascii_uppercase() == seq2[j - 1].to_ascii_uppercase() {
                aligned_identity += 1.0;
            }
            aligned_length += 1;
            i -= 1;
            j -= 1;
        } else if i > 0 && (score_matrix[i][j] - ins_matrix[i][j]).abs() < 1e-6 {
            aligned_seq1.push(seq1[i - 1]);
            aligned_seq2.push(b'-');
            i -= 1;
        } else if j > 0 && (score_matrix[i][j] - del_matrix[i][j]).abs() < 1e-6 {
            aligned_seq1.push(b'-');
            aligned_seq2.push(seq2[j - 1]);
            j -= 1;
        } else {
            break;
        }
    }

    aligned_seq1.reverse();
    aligned_seq2.reverse();

    let mut markup = String::with_capacity(aligned_seq1.len());
    for (&a, &b) in aligned_seq1.iter().zip(aligned_seq2.iter()) {
        let ch = if a.to_ascii_uppercase() == b.to_ascii_uppercase() && a != b'-' && b != b'-' {
            '|'
        } else if a == b'-' || b == b'-' {
            ' '
        } else if score_fn(a, b) > 0.0 {
            ':'
        } else {
            '.'
        };
        markup.push(ch);
    }

    AlignmentResult {
        aligned_seq1: String::from_utf8(aligned_seq1).unwrap(),
        aligned_seq2: String::from_utf8(aligned_seq2).unwrap(),
        aligned_length,
        aligned_identity: if aligned_length > 0 {
            aligned_identity / aligned_length as f64
        } else {
            0.0
        },
        score: score_matrix[len1][len2],
        alignment_markup: markup,
    }
}

pub fn needleman_wunsch_blosum62_internal(
    seq1: &str,
    seq2: &str,
    gap_open: f64,
    gap_extend: f64,
) -> AlignmentResult {
    needleman_wunsch_with_matrix(seq1, seq2, gap_open, gap_extend, |a, b| {
        blosum62_score(a, b) as f64
    })
}

pub fn needleman_wunsch_internal(
    seq1: &str,
    seq2: &str,
    match_score: f64,
    mismatch_penalty: f64,
    gap_open: f64,
    gap_extend: f64,
) -> AlignmentResult {
    needleman_wunsch_with_matrix(seq1, seq2, gap_open, gap_extend, |a, b| {
        if a.to_ascii_uppercase() == b.to_ascii_uppercase() {
            match_score
        } else {
            mismatch_penalty
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn align_identical_sequences() {
        let r = smith_waterman_internal("GATTACA", "GATTACA", 2.0, -1.0, -1.0, -0.5);
        assert_eq!(r.aligned_seq1, "GATTACA");
        assert_eq!(r.aligned_seq2, "GATTACA");
        assert_eq!(r.aligned_length, 7);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
        assert_eq!(r.alignment_markup, "|||||||");
    }

    #[test]
    fn align_acacacta_agcacaca() {
        let r = smith_waterman_internal("ACACACTA", "AGCACACA", 2.0, -1.0, -1.0, -0.5);
        assert_eq!(r.aligned_seq1, "A-CACACTA");
        assert_eq!(r.aligned_seq2, "AGCACAC-A");
        assert_eq!(r.aligned_length, 7);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn align_gattaca_gcatgcu() {
        let r = smith_waterman_internal("GATTACA", "GCATGCU", 2.0, -1.0, -1.0, -0.5);
        assert_eq!(r.aligned_seq1, "G-AT---TACA");
        assert_eq!(r.aligned_seq2, "GCATGCU----");
        assert_eq!(r.aligned_length, 3);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn blosum62_identical_sequences() {
        let r = smith_waterman_blosum62_internal("GATTACA", "GATTACA", -1.0, -0.5);
        assert_eq!(r.aligned_seq1, "GATTACA");
        assert_eq!(r.aligned_seq2, "GATTACA");
        assert_eq!(r.aligned_length, 7);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn smith_waterman_case_insensitive() {
        let r = smith_waterman_internal("gAttAcA", "GATTACA", 2.0, -1.0, -1.0, -0.5);
        assert_eq!(r.aligned_seq1, "gAttAcA");
        assert_eq!(r.aligned_seq2, "GATTACA");
        assert_eq!(r.aligned_length, 7);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
        assert_eq!(r.alignment_markup, "|||||||");
    }

    #[test]
    fn nw_identical_sequences() {
        let r = needleman_wunsch_internal("GATTACA", "GATTACA", 2.0, -1.0, -1.0, -0.5);
        assert_eq!(r.aligned_seq1, "GATTACA");
        assert_eq!(r.aligned_seq2, "GATTACA");
        assert_eq!(r.aligned_length, 7);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn needleman_wunsch_case_insensitive() {
        let r = needleman_wunsch_internal("GATTACA", "gattaca", 2.0, -1.0, -1.0, -0.5);
        assert_eq!(r.aligned_seq1, "GATTACA");
        assert_eq!(r.aligned_seq2, "gattaca");
        assert_eq!(r.aligned_length, 7);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
        assert_eq!(r.alignment_markup, "|||||||");
    }

    #[test]
    fn nw_blosum62_identical_sequences() {
        let r = needleman_wunsch_blosum62_internal("GATTACA", "GATTACA", -1.0, -0.5);
        assert_eq!(r.aligned_seq1, "GATTACA");
        assert_eq!(r.aligned_seq2, "GATTACA");
        assert_eq!(r.aligned_length, 7);
        assert!((r.aligned_identity - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn blosum62_symmetry() {
        for i in 0..20 {
            for j in 0..20 {
                assert_eq!(
                    BLOSUM62_MATRIX[i][j],
                    BLOSUM62_MATRIX[j][i],
                    "BLOSUM62 matrix is not symmetric at indices ({}, {})",
                    i, j
                );
            }
        }
    }

    #[test]
    fn sw_known_cases() {
        let seq1 = "MTFSSTSSAPPPSPLLPATRITVYGCGRDEAALFRRTAPRFGVEATLTEAAVSEENAEMAAGNQCISIDHKTPVTPATLRALHRAGVTYISTRSIGYNHIDVTYAAGVGISVENVTYSPAGVADYTLMLMLMAVRNAKSTVRRAELHDYRLNEIRGKELRDLTVGVIGTGRIGAAVVDRLRGFGSRVLAYGKRPTIAADYVSLDELLRSSDIVSLHVPLTPDTHHLLDQSRIRRMKSGAFVINTGRGPLIDTEALVPALESGRLSGAALDVIEGEEGIFYADCRNRTIESTWLPRLQKMPNVLISPHTAYYTDHALMDTVENSIINCLNFGSRKQHGVGQVGQVEGRHRIRGLFRRTRRFRQVRPGGRTQPRHREVPAVLRGDHEGRRLETLRRARPGLGERRLPS";
        let seq2 = "MSYRDLGLIDSEVIAERRVRALDDSSPSAVPTTGVRVFGCGHDEAVLFREMGTRLGITPSITEEAISETNAELARGNRCISVSHKTQIDNSTLLALSRVGVEYISTRSVGYNHIDVEFAASIGISVGNVDYSPDSVGDYTLMLMLMTVRHAKSIVRRADTHDYRLNDTRGRELRDLTVGVIGTGRIGTAVIDRLQGFGCRVLAHDSGPHASADYVPLDELLRQSDIVTLHTPLTADTHHLLDRQRIDQMKHGAYIVNTGRGPLLDTEALLSALESGRLGGAALDVVEGEEGIFYADCRNRLIENKALVRLQRLPNVLISPHSAYYTDHALNDTVENSLVNCLNFESGRTA";

        // tested againt EMBL-EBI alignment tool & against Biopython
        assert_eq!(
            smith_waterman_blosum62_internal(seq1, seq2, -10.0, -0.5).score,
            1178.0);

        let r5_005 = smith_waterman_blosum62_internal(seq1, seq2, -5.0, -0.05);
        assert!((r5_005.score - 1198.1).abs() < 0.01, "Expected score close to 1198.1, got {}", r5_005.score);

        assert_eq!(
            smith_waterman_blosum62_internal(seq1, seq2, -5.0, -0.5).score,
            1185.0);

        // NW Tested againt Biopython's implementation
        assert_eq!(
            needleman_wunsch_blosum62_internal(seq1, seq2, -10.0, -0.5).score,
            1130.0);

        assert_eq!(
            needleman_wunsch_blosum62_internal(seq1, seq2, -2.0, -0.5).score,
            1181.0);

        assert_eq!(
            needleman_wunsch_blosum62_internal(&seq1[10..100], &seq2[10..100], -2.0, -0.5).score,
            201.0);

        assert_eq!(
            needleman_wunsch_blosum62_internal(&seq1[10..100], &seq2[10..100], -7.0, -0.5).score,
            178.5);

    }
}
