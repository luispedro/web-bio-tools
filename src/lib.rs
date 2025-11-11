use serde_wasm_bindgen::to_value;
use wasm_bindgen::prelude::*;
use wasm_bindgen::JsValue;

mod alignment;
mod hmm;
#[cfg(all(feature = "python", not(target_arch = "wasm32")))]
mod python;

pub use alignment::AlignmentResult;
pub use hmm::HmmProfile;

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

#[wasm_bindgen]
pub fn parse_hmm(text: &str) -> Result<JsValue, JsValue> {
    let profile = hmm::parse_hmm(text).map_err(|e| JsValue::from_str(&e))?;
    to_value(&profile)
        .map_err(|err| JsValue::from_str(&format!("Failed to serialise HMM profile: {err}")))
}
