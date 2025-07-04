# Web Bio Tools

Rust implementation of some simple bioinformatics tools designed to run in the browser (using WebAssembly).

Live version: [https://web-bio-tools.netlify.app/](https://web-bio-tools.netlify.app/)


## Local development

For building, you can use `wasm-pack` to compile the Rust code to WebAssembly.
Make sure you have it installed:

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

