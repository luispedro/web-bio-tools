#!/usr/bin/env bash

set -ve
if test -d dist; then
  rm -rf dist.bak
  mv dist dist.bak
fi

wasm-pack build --target web
mkdir -p dist
cp -r pkg/ dist/
cp -r static/ dist/
cp -r index.html dist/
cp -r sidebar.html matrix.html hmm.html faq.html fna2faa.html deltavis.html deltavis.js changelog.html dist/
cp -r demo-data/ dist/
