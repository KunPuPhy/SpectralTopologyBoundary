#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/.."
rm -f Note.aux Note.log Note.out Note.toc Note.bbl Note.blg \
      Note.fls Note.fdb_latexmk Note.synctex.gz "Note.synctex(busy)" \
      NoteNotes.bib missfont.log
