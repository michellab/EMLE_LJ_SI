#!/bin/bash
# Build postg with postg.patch (external mesh and atomic weights). Usage: ./build_postg.sh <build directory>
set -euo pipefail
PATCH=$(readlink -f "$(dirname "$0")/postg.patch")
git clone https://github.com/aoterodelaroza/postg.git "$1"
cd "$1"
git checkout 1c227da0fa5490e44e3b2a5b0c9156d5aaf22cc5
git apply "$PATCH"
make
