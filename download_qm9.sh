#!/usr/bin/env bash
# Download QM9 (gdb9.sdf, gdb9.sdf.csv, uncharacterized.txt) into data/qm9/, or into $QM9_DIR if set.
set -euo pipefail
DIR="${QM9_DIR:-$(dirname "$0")/data/qm9}"
mkdir -p "$DIR"
curl -L https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/gdb9.tar.gz | tar -xz -C "$DIR"
curl -L -o "$DIR/uncharacterized.txt" https://ndownloader.figshare.com/files/3195404
ls -la "$DIR"
