#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="$SCRIPT_DIR"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$ROOT_DIR"

xelatex -interaction=nonstopmode -output-directory="$OUT_DIR" main.tex
cd "$OUT_DIR"
cp "$ROOT_DIR/ref.bib" . 2>/dev/null || true
bibtex main > /dev/null
cd "$ROOT_DIR"
xelatex -interaction=nonstopmode -output-directory="$OUT_DIR" main.tex
xelatex -interaction=nonstopmode -output-directory="$OUT_DIR" main.tex
cp "$OUT_DIR/main.pdf" "$ROOT_DIR/main.pdf"