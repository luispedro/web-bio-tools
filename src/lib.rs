use wasm_bindgen::prelude::*;
use serde_wasm_bindgen::to_value;
use wasm_bindgen::JsValue;

mod alignment;

pub use alignment::AlignmentResult;

#[wasm_bindgen]
pub fn smith_waterman(seq1: &str, seq2: &str) -> JsValue {
    smith_waterman_custom(seq1, seq2, 2, -1, -1)
}

#[wasm_bindgen]
pub fn smith_waterman_custom(
    seq1: &str,
    seq2: &str,
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> JsValue {
    let result = alignment::smith_waterman_internal(
        seq1,
        seq2,
        match_score,
        mismatch_penalty,
        gap_penalty,
    );
    to_value(&result).unwrap()
}

#[wasm_bindgen]
pub fn smith_waterman_blosum62(seq1: &str, seq2: &str, gap_penalty: i32) -> JsValue {
    let result = alignment::smith_waterman_blosum62_internal(seq1, seq2, gap_penalty);
    to_value(&result).unwrap()
}

#[wasm_bindgen]
pub fn needleman_wunsch(seq1: &str, seq2: &str) -> JsValue {
    needleman_wunsch_custom(seq1, seq2, 2, -1, -1)
}

#[wasm_bindgen]
pub fn needleman_wunsch_custom(
    seq1: &str,
    seq2: &str,
    match_score: i32,
    mismatch_penalty: i32,
    gap_penalty: i32,
) -> JsValue {
    let result = alignment::needleman_wunsch_internal(
        seq1,
        seq2,
        match_score,
        mismatch_penalty,
        gap_penalty,
    );
    to_value(&result).unwrap()
}

#[wasm_bindgen]
pub fn needleman_wunsch_blosum62(seq1: &str, seq2: &str, gap_penalty: i32) -> JsValue {
    let result = alignment::needleman_wunsch_blosum62_internal(seq1, seq2, gap_penalty);
    to_value(&result).unwrap()
}
