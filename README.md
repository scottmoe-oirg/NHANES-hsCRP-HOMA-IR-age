# Insulin Resistance and Systemic Inflammation Below the HbA1c Prediabetes Threshold

## Overview

This repository contains the manuscript, supplementary material, reproducible
analysis pipeline, computational outputs, and independent cross-software
validation for:

**Insulin Resistance and Systemic Inflammation Below the HbA1c Prediabetes
Threshold: Age-Dependent Associations in U.S. Adults, NHANES 2015–March 2020**

The study examines the association between insulin resistance, estimated using
the Homeostatic Model Assessment of Insulin Resistance (HOMA-IR), and systemic
inflammation, measured using high-sensitivity C-reactive protein (hs-CRP), among
U.S. adults without diagnosed diabetes and with HbA1c <5.7%.

The primary analysis tests whether the adjusted HOMA-IR–hs-CRP association
varies across adulthood.

## Study design

The analysis uses nationally representative data from:

- NHANES 2015–2016
- NHANES 2017–March 2020 pre-pandemic

The final analytic sample contains 3,079 adults aged 20 years or older with:

- no diagnosed diabetes;
- HbA1c <5.7%;
- a qualifying fasting-subsample weight;
- valid fasting glucose and insulin measurements;
- positive HOMA-IR and hs-CRP values; and
- complete primary covariate data.

The primary outcome is log-transformed hs-CRP. The primary exposure is
HOMA-IR, calculated as:

    HOMA-IR = fasting insulin × fasting glucose / 405

The primary survey-weighted regression model evaluates the interaction between
HOMA-IR and age group (20–39, 40–54, and ≥55 years), adjusting for HbA1c,
waist circumference, continuous age, sex, and non-HDL cholesterol, with
period-specific nuisance relationships.

The primary age interaction was:

    F(2,40) = 5.686345
    p = 0.006707391

Age-specific adjusted HOMA-IR coefficients were:

    20–39 years:  β =  0.026629
    40–54 years:  β =  0.051127
    ≥55 years:    β = -0.018904

The manuscript discusses the interpretation and limitations of these estimates,
including sensitivity to model specification and the cross-sectional design.

## Repository structure

### `hscrp_homair_reproducible_analysis/`

Canonical Python implementation of the analysis.

The pipeline reconstructs the analytic sample from the original NHANES
public-use XPT files, applies the survey design and combined-period fasting
weights, fits the primary model, produces age-specific estimates, and generates
prespecified validation and audit outputs.

Run from this directory with:

    python3 run_analysis.py

See the README within that directory for detailed instructions.

### `validation/`

Independent cross-software implementation in R using the `survey` package.

Canonical validation script:

    validate_analysis_R_v5.R

From the repository root, run:

    ./validation/run_validation.sh

The validation independently reconstructs the analytic sample from the raw
NHANES files and reproduces the primary analysis and selected sensitivity
analyses.

Reference validation output is provided in:

    R_validate_full_v5.text

### `manuscript/`

Contains the canonical manuscript, supplementary material, figures, tables,
aggregate statistical outputs, and supporting analyses.

Canonical manuscript:

    NHANES_hsCRP_HOMAIR_age_manuscript.tex
    NHANES_hsCRP_HOMAIR_age_manuscript.pdf

Canonical supplement:

    NHANES_hsCRP_HOMAIR_age_supplement.tex
    NHANES_hsCRP_HOMAIR_age_supplement.pdf

## Data

All source data used in this study are publicly available through the National
Center for Health Statistics (NCHS) NHANES program.

Raw NHANES XPT files are intentionally not included in this repository.
Derived participant-level analytic and quality-control datasets are also not
redistributed.

To reproduce the analysis, obtain the required public-use NHANES files and
place them in:

    data_raw/

The analysis pipeline then reconstructs the analytic dataset from those source
files.

## Validated Python environment

The reproducible analysis was successfully validated in a clean Python virtual
environment using:

- Python 3.13.9
- NumPy 2.3.4
- pandas 3.0.5
- patsy 1.0.2
- SciPy 1.18.0
- statsmodels 0.14.6

Exact package versions are recorded in:

    hscrp_homair_reproducible_analysis/requirements.txt

A clean-environment test beginning with the documented dependencies and raw
NHANES XPT files reproduced the locked analytic sample, survey-design degrees
of freedom, primary interaction test, and age-specific HOMA-IR coefficients,
with the pipeline reporting:

    Validation: PASS

These versions document the environment in which the analysis was validated;
they are not intended to imply that other software versions are incompatible.

## Independent R validation

The primary analysis was additionally implemented independently in R as a
cross-software validation.

The validated environment used:

- R 4.6.1
- survey 4.5
- haven 2.5.5

The independent implementation reproduced the primary analysis, crude model,
and log-HOMA-IR sensitivity analysis from the raw NHANES files.

## Reproducibility

The Python implementation is the canonical computational analysis.

The R implementation was developed as an independent cross-software validation
rather than as a second canonical pipeline. Agreement between the two
implementations provides an additional check on analytic-sample construction,
NHANES survey-design handling, model specification, and statistical inference.

Raw source data and participant-level derived datasets are intentionally
excluded from version control. Aggregate statistical results and computational
audit outputs required to document the analyses are retained.

## Authors

**Scott Moerschbacher** — Open Inquiry Research Group

**Janice Baker** — Open Inquiry Research Group; Department of Sociology, University of Nevada, Reno

## Citation

A versioned archival release and DOI will be provided through Zenodo.

Citation metadata will also be provided in `CITATION.cff`.

## License

The computational code in this repository is distributed under the MIT License; see `LICENSE`.

The NHANES source data are distributed separately by the National Center for
Health Statistics and remain subject to their applicable terms and
documentation.
