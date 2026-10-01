#!/usr/bin/env bash
# Build the Agent Master fixture-execution image and print its immutable ID.
#
# Only git-tracked package files at HEAD enter the build context, so ignored
# files (credentials such as .env.local, local data, outputs) cannot leak into
# the image. Set PYCSAMT_FIXTURE_IMAGE to the printed sha256 ID.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
git_() { git -c safe.directory='*' -C "$here" "$@"; }
repo=$(git_ rev-parse --show-toplevel)
rev=$(git_ rev-parse --short HEAD)

context=$(mktemp -d)
trap 'rm -rf "$context"' EXIT
cp "$here/Dockerfile" "$context/"
mkdir "$context/src"
git -c safe.directory='*' -C "$repo" archive HEAD \
    pycsamt pyproject.toml README.md LICENSE.md MANIFEST.in \
    | tar -x -C "$context/src"

docker build --label org.pycsamt.revision="$rev" \
    -t "pycsamt-agent-fixture:$rev" "$context" >&2
docker image inspect --format '{{.Id}}' "pycsamt-agent-fixture:$rev"
