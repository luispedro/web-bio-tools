use wasm_bindgen::prelude::*;
use serde_wasm_bindgen::to_value;
use wasm_bindgen::JsValue;

mod alignment;

pub use alignment::AlignmentResult;

#[wasm_bindgen]
pub fn smith_waterman(seq1: &str, seq2: &str) -> JsValue {
    let result = alignment::smith_waterman_internal(seq1, seq2);
    to_value(&result).unwrap()
}
