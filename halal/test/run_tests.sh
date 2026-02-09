#!/bin/bash
set -e

cd "$(dirname "$0")/.."

echo "Building HalalSeq..."
make clean
make all
make test
