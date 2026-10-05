#!/bin/bash
set -euo pipefail

# ============================================================
# Independent R cross-software validation
# hsCRP–HOMA-IR NHANES study
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DATA_DIR="${PROJECT_ROOT}/data_raw"
R_SCRIPT="${SCRIPT_DIR}/validate_analysis_R_v5.R"
OUTPUT_FILE="${SCRIPT_DIR}/R_validate_full_v5.text"

echo "============================================================"
echo " hsCRP–HOMA-IR Independent R Validation"
echo "============================================================"
echo
echo "Validation script:"
echo "  ${R_SCRIPT}"
echo
echo "Raw data directory:"
echo "  ${DATA_DIR}"
echo
echo "Output:"
echo "  ${OUTPUT_FILE}"
echo

if [ ! -f "${R_SCRIPT}" ]; then
    echo "ERROR: R validation script not found."
    exit 1
fi

if [ ! -d "${DATA_DIR}" ]; then
    echo "ERROR: Raw data directory not found."
    exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript was not found on PATH."
    exit 1
fi

echo "Running validation..."
echo

Rscript "${R_SCRIPT}" \
    --data-dir "${DATA_DIR}" \
    > "${OUTPUT_FILE}"

echo "Validation completed successfully."
echo
echo "Results written to:"
echo "  ${OUTPUT_FILE}"
echo
echo "============================================================"
