# Web Bio Tools
[![Tests](https://github.com/luispedro/web-bio-tools/actions/workflows/test.yml/badge.svg?branch=main)](https://github.com/luispedro/web-bio-tools/actions/workflows/test.yml) [![Deployed](https://img.shields.io/website?url=https%3A%2F%2Fweb-bio-tools.big-data-biology.org&label=deployed)](https://web-bio-tools.big-data-biology.org/)

**Warning: This project is experimental and untested.**

Rust implementation of some simple bioinformatics tools designed to run in the browser (using WebAssembly).

Live version: [https://web-bio-tools.big-data-biology.org/](https://web-bio-tools.big-data-biology.org/)


## Local development

For building, you can use `wasm-pack` to compile the Rust code to WebAssembly.
Make sure you have it installed (you can install it via `cargo install
wasm-pack`) and then run:

```bash
wasm-pack build --target web
```

You can test the tools in your browser by running a local server. For example, using Python:

```bash
python -m http.server
```

Then, open your browser and navigate to
[http://localhost:8000](http://localhost:8000).


## License

This project is licensed under the MIT License. See the
[LICENSE.MIT](LICENSE.MIT) file for details.

## Authors

- [Luis Pedro Coelho](https://luispedro.org/)

