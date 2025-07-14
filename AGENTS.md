# Repo AGENTS Instructions

This project contains Rust code compiled to WebAssembly and a small Python package.

## Formatting
- Run `cargo fmt` on all Rust files before committing.

## Testing
- Run `cargo test` and ensure it succeeds before creating a PR.
- If you modify Python code or tests, build the Python extension with `maturin develop` (or `pip install -e .`) and run `pytest`.

## Pull Request
- Summarize the main changes and include the results of `cargo test` (and `pytest` if run) in the PR description.
