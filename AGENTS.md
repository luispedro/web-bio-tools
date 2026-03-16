# Repo AGENTS Instructions

This project contains Rust code compiled to WebAssembly and a small Python package.

## Project Overview

Web Bio Tools is a Rust bioinformatics toolkit compiled to WebAssembly for browser execution, with optional Python bindings (PyO3). It provides three tools: sequence alignment (Smith-Waterman/Needleman-Wunsch), DNA-to-protein translation (6-frame), and HMM profile viewer (HMMER3 format with forward-backward algorithm). Runs entirely client-side — no backend server.

Live: https://web-bio-tools.big-data-biology.org/

## Build Commands

```bash
# WebAssembly build (required for browser)
wasm-pack build --target web

# Python bindings (optional, requires 'python' feature)
maturin develop          # development install
maturin build --features python  # release build

# Local dev server
python -m http.server
```

## Testing

```bash
# Rust tests (must pass before PRs)
cargo test --verbose

# Python tests (if Python code changed)
maturin develop && pytest
```

## Formatting

Run `cargo fmt` on all Rust files before committing.

## Architecture

All core algorithms are in `src/`:
- **`lib.rs`** — WASM entry points (`#[wasm_bindgen]` exports), ties modules together
- **`alignment.rs`** — Smith-Waterman (local) and Needleman-Wunsch (global) alignment with BLOSUM62 or custom scoring; affine gap penalties
- **`hmm.rs`** — HMMER3 HMM parsing and forward-backward probability computation
- **`fna2faa.rs`** — Codon table and DNA→protein translation (handles IUPAC ambiguous nucleotides via bitmask)
- **`translation.rs`** — Frame translation interface (all 6 reading frames)
- **`python.rs`** — PyO3 bindings (behind `python` feature flag)

Frontend is plain HTML/JS (no framework): `index.html` (alignment), `fna2faa.html` (translator), `hmm.html` (HMM viewer). They load WASM from `pkg/` and use Bootstrap 4 + jQuery.

## Key Dependencies

- `wasm-bindgen` / `serde-wasm-bindgen` — Rust↔JS interop
- `pyo3` — Python FFI (optional feature `python`)
- Python tests use `pytest`, `hypothesis` (property-based), and `biopython` as reference

## CI

GitHub Actions (`.github/workflows/test.yml`): runs `cargo test --verbose` on push to main and all PRs. Python tests are not in CI.

## PR Requirements

Summarize changes and include `cargo test` output (and `pytest` if applicable) in PR description.
