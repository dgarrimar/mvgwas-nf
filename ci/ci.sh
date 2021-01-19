#!/bin/bash
set -e
set -o pipefail

OUT_DIR=ci

case "$1" in
  run)
    shift
    echo "Running test pipeline..." >&2
    nextflow run . -resume -with-docker -w ${OUT_DIR}/work --dir ${OUT_DIR} $@
    ;;
  validate)
    echo "Validating test results..." >&2
    spiff ${OUT_DIR}/mvgwas.tsv ${OUT_DIR}/mvgwas.ref.tsv
    ;;
  cleanup)
    echo "Cleaning up test results..." >&2
    rm -rf ${OUT_DIR}/mvgwas.tsv ${OUT_DIR}/work
    ;;
  *)
    echo "Usage: ci.sh {run|validate|cleanup}" >&2
    exit 1
esac
