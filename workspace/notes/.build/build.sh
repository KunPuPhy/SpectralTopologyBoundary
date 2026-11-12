#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="$SCRIPT_DIR"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
cd "$ROOT_DIR"

xelatex -interaction=nonstopmode -output-directory="$OUT_DIR" Note.tex
cd "$OUT_DIR"
bibtex Note > /dev/null
cd "$ROOT_DIR"
xelatex -interaction=nonstopmode -output-directory="$OUT_DIR" Note.tex
xelatex -interaction=nonstopmode -output-directory="$OUT_DIR" Note.tex
cp "$OUT_DIR/Note.pdf" "$ROOT_DIR/Note.pdf"
