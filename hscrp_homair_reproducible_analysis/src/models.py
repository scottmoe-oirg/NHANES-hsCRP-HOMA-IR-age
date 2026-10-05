from __future__ import annotations

import pandas as pd
from .survey import survey_linear_regression, joint_wald_test, linear_contrast


def _primary_rhs(covariates, period_specific_nuisance=True):
    cov = " + ".join(covariates)
    if period_specific_nuisance:
        return f"HOMA_IR * C(age_group) + C(period) * ({cov})"
    return f"HOMA_IR * C(age_group) + C(period) + {cov}"


def interaction_summary(result):
    interaction_terms = [n for n in result["names"] if n.startswith("HOMA_IR:C(age_group)")]
    test = joint_wald_test(result, interaction_terms)
    ref = "HOMA_IR"
    slopes = []
    for label, extra in [
        ("20-39", None),
        ("40-54", "HOMA_IR:C(age_group)[T.40-54]"),
        ("55+", "HOMA_IR:C(age_group)[T.55+]"),
    ]:
        coeff = {ref: 1.0}
        if extra is not None:
            coeff[extra] = 1.0
        row = linear_contrast(result, coeff)
        row["age_group"] = label
        row["n_model"] = result["n"]
        row["design_df"] = result["df"]
        slopes.append(row)
    return test, pd.DataFrame(slopes)


def fit_primary(data, domain, config):
    covariates = list(config["variables"]["primary_covariates"])
    rhs = _primary_rhs(covariates, bool(config["survey"]["period_specific_nuisance"]))
    formula = f"log_hs_CRP ~ {rhs}"
    result = survey_linear_regression(data, formula, domain)
    test, slopes = interaction_summary(result)
    return formula, result, test, slopes


def fit_crude(data, domain):
    # Survey period is retained only as a combined-period nuisance intercept.
    formula = "log_hs_CRP ~ HOMA_IR * C(age_group) + C(period)"
    result = survey_linear_regression(data, formula, domain)
    test, slopes = interaction_summary(result)
    return formula, result, test, slopes


def fit_with_covariates(data, domain, covariates, period_specific=True):
    if len(covariates) == 0:
        return fit_crude(data, domain)
    cov = " + ".join(covariates)
    if period_specific:
        formula = f"log_hs_CRP ~ HOMA_IR * C(age_group) + C(period) * ({cov})"
    else:
        formula = f"log_hs_CRP ~ HOMA_IR * C(age_group) + C(period) + {cov}"
    result = survey_linear_regression(data, formula, domain)
    test, slopes = interaction_summary(result)
    return formula, result, test, slopes
