#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
cd "$SCRIPT_DIR"

VERSION="3.1.3"
ARCHIVE_NAME="iqtree-${VERSION}-macOS.zip"
RELEASE_DIR="iqtree-${VERSION}-macOS"
URL="https://github.com/iqtree/iqtree3/releases/download/v${VERSION}/${ARCHIVE_NAME}"
EXPECTED_SHA256="3ede720b2871205cfc322580c9491cadd152bbedb3e1e9dd3e559951ed6d28c7"
TOOLS_DIR="tools"
ARCHIVE_PATH="${TOOLS_DIR}/${ARCHIVE_NAME}"
IQTREE_BIN="${TOOLS_DIR}/${RELEASE_DIR}/bin/iqtree3"

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "This installer is for macOS only." >&2
  exit 2
fi

mkdir -p "$TOOLS_DIR" tmp

if [[ ! -f "$ARCHIVE_PATH" ]]; then
  echo "Downloading IQ-TREE ${VERSION} macOS Universal from the official release..."
  curl -L --fail --show-error "$URL" -o "$ARCHIVE_PATH"
fi

printf '%s  %s\n' "$EXPECTED_SHA256" "$ARCHIVE_PATH" | shasum -a 256 -c -

if [[ ! -x "$IQTREE_BIN" ]]; then
  unzip -q -o "$ARCHIVE_PATH" -d "$TOOLS_DIR"
  chmod +x "$IQTREE_BIN"
fi

VERSION_TEXT="$($IQTREE_BIN --version 2>&1 | head -n 1)"
echo "$VERSION_TEXT"
if [[ "$VERSION_TEXT" == *"single-core"* ]]; then
  echo "The downloaded binary unexpectedly reports a single-core build." >&2
  exit 3
fi

TEST_PREFIX="tmp/iqtree_multicore_install_test"
"$IQTREE_BIN" \
  -s "${TOOLS_DIR}/${RELEASE_DIR}/example.phy" \
  -m JC \
  -nt 2 \
  -pre "$TEST_PREFIX" \
  -redo >/dev/null

if ! grep -Eq '2 threads|2 CPU cores' "${TEST_PREFIX}.log"; then
  echo "IQ-TREE ran, but the two-thread verification was not found in its log." >&2
  exit 4
fi

echo "Multicore verification passed: 2 threads."
echo "Installed executable: ${SCRIPT_DIR}/${IQTREE_BIN}"
echo "phylogeny_workflow.py will select this project-local executable automatically."
