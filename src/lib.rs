use serde_wasm_bindgen::to_value;
use wasm_bindgen::prelude::*;
use wasm_bindgen::JsValue;

mod alignment;
mod fna2faa;
mod hmm;
#[cfg(all(feature = "python", not(target_arch = "wasm32")))]
mod python;
mod translation;

pub use alignment::AlignmentResult;
pub use translation::{translate_all_frames, translate_frame};

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
    let hmm = hmm::parse_hmm(text).map_err(|err| JsValue::from_str(&err))?;
    to_value(&hmm).map_err(|err| JsValue::from_str(&format!("Failed to serialize HMM: {}", err)))
}

#[wasm_bindgen]
pub fn translate_dna_frame(
    seq: &str,
    frame: usize,
    stop_at_first_stop: bool,
) -> Result<String, JsValue> {
    translation::translate_frame(seq, frame, stop_at_first_stop)
        .map_err(|err| JsValue::from_str(&err))
}

#[wasm_bindgen]
pub fn translate_dna_all_frames(seq: &str, stop_at_first_stop: bool) -> Result<JsValue, JsValue> {
    let translations = translation::translate_all_frames(seq, stop_at_first_stop)
        .map_err(|err| JsValue::from_str(&err))?;
    to_value(&translations)
        .map_err(|err| JsValue::from_str(&format!("Failed to serialize translations: {}", err)))
}
