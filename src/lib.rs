use serde_wasm_bindgen::to_value;
use wasm_bindgen::prelude::*;
use wasm_bindgen::JsValue;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyDict;
#[cfg(feature = "python")]
use pyo3::{wrap_pyfunction, PyObject};

mod alignment;

pub use alignment::AlignmentResult;

#[wasm_bindgen]
pub fn smith_waterman(seq1: &str, seq2: &str) -> JsValue {
    smith_waterman_custom(seq1, seq2, 2.0, -1.0, -1.0, -0.5)
}

#[wasm_bindgen]
pub fn smith_waterman_custom(
    seq1: &str,
    seq2: &str,
    match_score: f64,
    mismatch_penalty: f64,
    gap_open: f64,
    gap_extend: f64,
) -> JsValue {
    let result = alignment::smith_waterman_internal(
        seq1,
        seq2,
        match_score,
        mismatch_penalty,
        gap_open,
        gap_extend,
    );
    to_value(&result).unwrap()
}

#[wasm_bindgen]
pub fn smith_waterman_blosum62(seq1: &str, seq2: &str, gap_open: f64, gap_extend: f64) -> JsValue {
    let result = alignment::smith_waterman_blosum62_internal(seq1, seq2, gap_open, gap_extend);
    to_value(&result).unwrap()
}

#[wasm_bindgen]
pub fn needleman_wunsch(seq1: &str, seq2: &str) -> JsValue {
    needleman_wunsch_custom(seq1, seq2, 2.0, -1.0, -1.0, -0.5)
}

#[wasm_bindgen]
pub fn needleman_wunsch_custom(
    seq1: &str,
    seq2: &str,
    match_score: f64,
    mismatch_penalty: f64,
    gap_open: f64,
    gap_extend: f64,
) -> JsValue {
    let result = alignment::needleman_wunsch_internal(
        seq1,
        seq2,
        match_score,
        mismatch_penalty,
        gap_open,
        gap_extend,
    );
    to_value(&result).unwrap()
}

#[wasm_bindgen]
pub fn needleman_wunsch_blosum62(
    seq1: &str,
    seq2: &str,
    gap_open: f64,
    gap_extend: f64,
) -> JsValue {
    let result = alignment::needleman_wunsch_blosum62_internal(seq1, seq2, gap_open, gap_extend);
    to_value(&result).unwrap()
}

#[cfg(feature = "python")]
fn result_to_dict(py: Python<'_>, r: AlignmentResult) -> PyObject {
    let d = PyDict::new(py);
    d.set_item("aligned_seq1", r.aligned_seq1).unwrap();
    d.set_item("aligned_seq2", r.aligned_seq2).unwrap();
    d.set_item("aligned_length", r.aligned_length).unwrap();
    d.set_item("aligned_identity", r.aligned_identity).unwrap();
    d.set_item("score", r.score).unwrap();
    d.set_item("alignment_markup", r.alignment_markup).unwrap();
    d.into()
}

// Python bindings
#[cfg(feature = "python")]
#[pyfunction]
fn py_smith_waterman(seq1: &str, seq2: &str, match_score: Option<f64>, mismatch_penalty: Option<f64>, gap_open: Option<f64>, gap_extend: Option<f64>) -> PyObject {
    Python::with_gil(|py| {
        let r = alignment::smith_waterman_internal(
            seq1,
            seq2,
            match_score.unwrap_or(2.0),
            mismatch_penalty.unwrap_or(-1.0),
            gap_open.unwrap_or(-1.0),
            gap_extend.unwrap_or(-0.5),
        );
        result_to_dict(py, r)
    })
}

#[cfg(feature = "python")]
#[pyfunction]
fn py_smith_waterman_blosum62(seq1: &str, seq2: &str, gap_open: Option<f64>, gap_extend: Option<f64>) -> PyObject {
    Python::with_gil(|py| {
        let r = alignment::smith_waterman_blosum62_internal(
            seq1,
            seq2,
            gap_open.unwrap_or(-1.0),
            gap_extend.unwrap_or(-0.5),
        );
        result_to_dict(py, r)
    })
}

#[cfg(feature = "python")]
#[pyfunction]
fn py_needleman_wunsch(seq1: &str, seq2: &str, match_score: Option<f64>, mismatch_penalty: Option<f64>, gap_open: Option<f64>, gap_extend: Option<f64>) -> PyObject {
    Python::with_gil(|py| {
        let r = alignment::needleman_wunsch_internal(
            seq1,
            seq2,
            match_score.unwrap_or(2.0),
            mismatch_penalty.unwrap_or(-1.0),
            gap_open.unwrap_or(-1.0),
            gap_extend.unwrap_or(-0.5),
        );
        result_to_dict(py, r)
    })
}

#[cfg(feature = "python")]
#[pyfunction]
fn py_needleman_wunsch_blosum62(seq1: &str, seq2: &str, gap_open: Option<f64>, gap_extend: Option<f64>) -> PyObject {
    Python::with_gil(|py| {
        let r = alignment::needleman_wunsch_blosum62_internal(
            seq1,
            seq2,
            gap_open.unwrap_or(-1.0),
            gap_extend.unwrap_or(-0.5),
        );
        result_to_dict(py, r)
    })
}

#[cfg(feature = "python")]
#[pymodule]
fn web_bio_tools(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_smith_waterman, m)?)?;
    m.add_function(wrap_pyfunction!(py_smith_waterman_blosum62, m)?)?;
    m.add_function(wrap_pyfunction!(py_needleman_wunsch, m)?)?;
    m.add_function(wrap_pyfunction!(py_needleman_wunsch_blosum62, m)?)?;
    Ok(())
}
