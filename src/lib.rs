use serde_wasm_bindgen::to_value;
use wasm_bindgen::prelude::*;
use wasm_bindgen::JsValue;

mod alignment;
mod fna2faa;
mod hmm;
#[cfg(all(feature = "python", not(target_arch = "wasm32")))]
mod python;

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

#[wasm_bindgen]
pub fn parse_hmm(text: &str) -> Result<JsValue, JsValue> {
    let hmm = hmm::parse_hmm(text).map_err(|err| JsValue::from_str(&err))?;
    to_value(&hmm).map_err(|err| JsValue::from_str(&format!("Failed to serialize HMM: {}", err)))
}

#[wasm_bindgen]
pub fn translate_frame(seq: &str, frame: i32, stop_at_first: bool) -> Result<JsValue, JsValue> {
    if frame < -3 || frame > 2 {
        return Err(JsValue::from_str("Frame must be between -3 and 2."));
    }

    let encoder = fna2faa::CodonEncoder::mk_encoder();
    let result = fna2faa::translate_frame_internal(&encoder, seq, frame as i8, stop_at_first)
        .map_err(|err| JsValue::from_str(&err))?;
    to_value(&result)
        .map_err(|err| JsValue::from_str(&format!("Failed to serialize translation: {}", err)))
}

#[wasm_bindgen]
pub fn translate_all_frames(seq: &str, stop_at_first: bool) -> Result<JsValue, JsValue> {
    let encoder = fna2faa::CodonEncoder::mk_encoder();
    let result = fna2faa::translate_all_frames_internal(&encoder, seq, stop_at_first);
    to_value(&result)
        .map_err(|err| JsValue::from_str(&format!("Failed to serialize translation: {}", err)))
}
