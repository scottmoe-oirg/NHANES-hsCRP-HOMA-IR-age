# Independent R cross-software validation — v5

This script independently reconstructs the locked NHANES analytic sample from
the raw public-use XPT files and validates three analyses with R's `survey`
package:

1. Primary adjusted model
2. Crude model
3. Adjusted log-HOMA-IR sensitivity model

For every joint interaction test, the script preserves and prints
`survey::regTermTest()`'s default result, then evaluates the same independently
obtained Wald F statistic using the prespecified NHANES survey design degrees
of freedom.

## Locked secondary targets

Crude model:
- F = 1.70699985
- p = 0.19435864
- slopes = 0.13079796, 0.15534031, 0.10172919

Log-HOMA sensitivity:
- F ≈ 1.86
- p ≈ 0.169

Only rounded values were retained for the log-HOMA sensitivity, so its automated
tolerances reflect that documented precision. The R age-specific log-HOMA
slopes are printed diagnostically rather than compared to invented extra digits.

## Validated software environment

The independent cross-software validation was successfully run using:

- R 4.6.1 (2026-06-24)
- survey 4.5
- haven 2.5.5

These versions document the software environment in which the validation
was confirmed to reproduce the locked Python analysis to numerical precision.
They should be interpreted as validated versions rather than strict minimum
software requirements.

## Run

```bash
Rscript validate_analysis_R_v5.R --data-dir /path/to/data_raw
```

A fully successful run ends with:

```text
FULL CROSS-SOFTWARE VALIDATION: PASS
Primary + crude + log-HOMA sensitivity reproduced.
```
