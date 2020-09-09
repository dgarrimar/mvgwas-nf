#!/bin/bash
set -e
set -o pipefail

OUT_DIR=ci

case "$1" in
  run)
    shift
    echo "Running test pipeline..." >&2
    nextflow run . -resume -with-docker --dir ${OUT_DIR} $@
    ;;
  validate)
    echo "Validating test results..." >&2
    md5sum -c ${OUT_DIR}/md5s.txt
    ;;
  cleanup)
    echo "Cleaning up test results..." >&2
    find ${OUT_DIR} -name mvgwas.tsv -exec rm {} \+
    ;;
  *)
    echo "Usage: ci.sh {run|validate}" >&2
    exit 1
esac
