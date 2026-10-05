# hsCRP-HOMA-IR NHANES reproducible analysis

This workspace is the canonical reproducibility backbone for the manuscript on the age-dependent association between HOMA-IR and hs-CRP below the HbA1c prediabetes threshold.

## Current reproducibility milestone

Starting from the public NHANES XPT files, the pipeline:

1. harmonizes NHANES 2015-2016 and 2017-March 2020 pre-pandemic data;
2. reconstructs the corrected analytic population (n = 3,079);
3. applies duration-adjusted fasting-subsample weights;
4. preserves the complex masked-stratum/PSU design and uses domain-analysis score zeroing;
5. reproduces the locked primary log(hs-CRP) HOMA-IR-by-age-group model;
6. validates the interaction and age-specific slopes against locked manuscript values;
7. reproduces the crude-vs-adjusted, single-covariate, drop-one, and sequential-adjustment audits.

## Raw XPT files required

Place the following files in `data_raw/`, or point `--data-dir` to the directory containing them:

2015-2016: `DEMO_I.xpt`, `BMX_I.xpt`, `DIQ_I.xpt`, `FASTQX_I.xpt`, `GHB_I.xpt`, `GLU_I.xpt`, `INS_I.xpt`, `HSCRP_I.xpt`, `TCHOL_I.xpt`, `HDL_I.xpt`

2017-March 2020: `P_DEMO.xpt`, `P_BMX.xpt`, `P_DIQ.xpt`, `P_FASTQX.xpt`, `P_GHB.xpt`, `P_GLU.xpt`, `P_INS.xpt`, `P_HSCRP.xpt`, `P_TCHOL.xpt`, `P_HDL.xpt`

The raw data are deliberately not bundled with this code package.

## Validated software environment

The reproducible analysis was successfully validated in a clean Python
virtual environment using:

- Python 3.13.9
- NumPy 2.3.4
- pandas 3.0.5
- patsy 1.0.2
- SciPy 1.18.0
- statsmodels 0.14.6

The exact package versions are recorded in `requirements.txt`.

A clean-environment test starting from the documented dependencies and raw
NHANES XPT files reproduced the locked analytic sample, survey design degrees
of freedom, primary interaction test, and age-specific HOMA-IR coefficients,
with the pipeline reporting `Validation: PASS`.

These versions document the environment in which the analysis was validated;
they are not intended to imply that other software versions are incompatible.

## Run

```bash
python3 run_analysis.py --data-dir /path/to/xpt/files --mode full_audit
```

Modes:

- `primary`: cohort construction + primary model + validation
- `crude_adjusted`: primary + crude-vs-adjusted audit
- `sequential_adjustment`: primary + single-covariate/drop-one/sequential audits
- `full_audit`: all currently implemented analyses

## Locked primary model

Outcome: natural log of hs-CRP.

Exposure/effect modification: `HOMA_IR * age_group`, with age groups 20-39, 40-54, and >=55 years.

Primary harmonized covariates: HbA1c, waist circumference, continuous age, sex, and non-HDL cholesterol.

For the combined-period primary model, survey period is permitted to have period-specific nuisance relationships with each primary adjustment covariate. HOMA-IR age-specific slopes remain the common primary estimands.

## Important fasting-weight correction

Some nominal zero fasting weights in the XPT files are encoded as tiny positive floating-point values (~5.4e-79). The configuration therefore defines a numerical zero tolerance. This prevents nominal zero-weight participants from being counted in the unweighted analytic n while leaving weighted estimates unchanged.

## Validation

`outputs/qc/validation.json` compares the run against the locked values:

- n = 3,079
- survey design df = 40
- HOMA-IR x age-group interaction F = 5.6863452144
- p = 0.00670739084
- slopes = 0.02662879194, 0.05112736866, -0.01890406625

The run raises an error if the locked values are not reproduced within the configured numerical tolerance.

## Scientific status of audits

The crude/adjusted and sequential/drop-one analyses are diagnostic/supportive analyses. They are kept separate from the locked primary inferential analysis so exploratory auditing cannot silently alter the primary model.
