use crate::AlignmentResult;
use pyo3::prelude::*;

#[pyclass]
pub struct PyAlignmentResult {
    #[pyo3(get)]
    pub aligned_seq1: String,
    #[pyo3(get)]
    pub aligned_seq2: String,
    #[pyo3(get)]
    pub aligned_length: usize,
    #[pyo3(get)]
    pub aligned_identity: f64,
    #[pyo3(get)]
    pub score: f64,
    #[pyo3(get)]
    pub alignment_markup: String,
}

impl From<AlignmentResult> for PyAlignmentResult {
    fn from(r: AlignmentResult) -> Self {
        Self {
            aligned_seq1: r.aligned_seq1,
            aligned_seq2: r.aligned_seq2,
            aligned_length: r.aligned_length,
            aligned_identity: r.aligned_identity,
            score: r.score,
            alignment_markup: r.alignment_markup,
        }
    }
}


#[pymodule]
mod web_bio_tools {
    use pyo3::prelude::*;
    use pyo3::exceptions::PyValueError;
    use crate::alignment;

    #[pymodule_export]
    use super::PyAlignmentResult;

    #[pyfunction]
    fn smith_waterman(seq1: &str, seq2: &str) -> PyAlignmentResult {
        alignment::smith_waterman_internal(seq1, seq2, 2.0, -1.0, -1.0, -0.5).into()
    }

    #[pyfunction]
    fn smith_waterman_custom(
        seq1: &str,
        seq2: &str,
        match_score: f64,
        mismatch_penalty: f64,
        gap_open: f64,
        gap_extend: f64,
    ) -> PyAlignmentResult {
        alignment::smith_waterman_internal(
            seq1,
            seq2,
            match_score,
            mismatch_penalty,
            gap_open,
            gap_extend,
        )
        .into()
    }

    #[pyfunction]
    fn smith_waterman_blosum62(
        seq1: &str,
        seq2: &str,
        gap_open: f64,
        gap_extend: f64,
    ) -> PyAlignmentResult {
        alignment::smith_waterman_blosum62_internal(seq1, seq2, gap_open, gap_extend).into()
    }

    #[pyfunction]
    fn needleman_wunsch(seq1: &str, seq2: &str) -> PyAlignmentResult {
        alignment::needleman_wunsch_internal(seq1, seq2, 2.0, -1.0, -1.0, -0.5).into()
    }

    #[pyfunction]
    fn needleman_wunsch_custom(
        seq1: &str,
        seq2: &str,
        match_score: f64,
        mismatch_penalty: f64,
        gap_open: f64,
        gap_extend: f64,
    ) -> PyAlignmentResult {
        alignment::needleman_wunsch_internal(
            seq1,
            seq2,
            match_score,
            mismatch_penalty,
            gap_open,
            gap_extend,
        )
        .into()
    }

    #[pyfunction]
    fn needleman_wunsch_blosum62(
        seq1: &str,
        seq2: &str,
        gap_open: f64,
        gap_extend: f64,
    ) -> PyAlignmentResult {
        alignment::needleman_wunsch_blosum62_internal(seq1, seq2, gap_open, gap_extend).into()
    }


    #[pyfunction]
    fn translate_dna_frame(seq: &str, frame: i8, stop_at_first_stop: bool) -> PyResult<String> {
        crate::translation::translate_frame(seq, frame, stop_at_first_stop)
            .map_err(|err| PyErr::new::<PyValueError, _>(err))
    }
}

