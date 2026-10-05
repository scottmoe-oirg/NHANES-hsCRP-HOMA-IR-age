from __future__ import annotations

import numpy as np
import pandas as pd
import patsy
from scipy import stats


def survey_linear_regression(data, formula, domain_mask,
                             weight="analysis_weight", strata="stratum_u", psu="psu_u"):
    """Survey-weighted linear regression with PSU/stratum Taylor-linearized covariance.

    Point estimates are weighted least squares. For variance estimation, score
    contributions are zeroed outside the analytic domain while the complete public-use
    survey design frame is retained. The design degrees of freedom are M-H, where M is
    the number of contributing PSUs and H the number of strata with >=2 PSUs.
    """
    y, X = patsy.dmatrices(formula, data, return_type="dataframe", NA_action="drop")
    idx = X.index
    dom = pd.Series(domain_mask, index=data.index).reindex(idx).fillna(False).astype(bool)
    good = (
        dom & data.loc[idx, weight].notna() & data.loc[idx, weight].gt(0)
        & data.loc[idx, strata].notna() & data.loc[idx, psu].notna()
    )
    Xd = X.loc[good].to_numpy(float)
    yd = y.loc[good].to_numpy(float).ravel()
    wd = data.loc[idx[good], weight].to_numpy(float)

    bread = Xd.T @ (wd[:, None] * Xd)
    bread_inv = np.linalg.inv(bread)
    beta = bread_inv @ (Xd.T @ (wd * yd))

    resid = yd - Xd @ beta
    p = Xd.shape[1]
    scores = pd.DataFrame(0.0, index=data.index, columns=range(p))
    scores.loc[idx[good], :] = (wd * resid)[:, None] * Xd

    tmp = data[[strata, psu]].copy()
    score_cols = []
    for j in range(p):
        col = f"s{j}"
        score_cols.append(col)
        tmp[col] = scores[j]
    psu_scores = tmp.groupby([strata, psu], dropna=True)[score_cols].sum().reset_index()

    meat = np.zeros((p, p))
    H = M = 0
    for _, g in psu_scores.groupby(strata):
        U = g[score_cols].to_numpy(float)
        m = len(U)
        if m < 2:
            continue
        H += 1
        M += m
        centered = U - U.mean(axis=0)
        meat += (m / (m - 1.0)) * centered.T @ centered

    cov = bread_inv @ meat @ bread_inv
    se = np.sqrt(np.diag(cov))
    design_df = M - H
    crit = stats.t.ppf(0.975, design_df)
    t = beta / se
    pval = 2 * stats.t.sf(np.abs(t), design_df)

    return {
        "names": list(X.columns), "beta": beta, "cov": cov, "se": se,
        "t": t, "p": pval, "ci_lo": beta - crit * se, "ci_hi": beta + crit * se,
        "df": int(design_df), "n": int(good.sum()),
    }


def joint_wald_test(result, terms):
    ids = [result["names"].index(t) for t in terms]
    b = result["beta"][ids]
    V = result["cov"][np.ix_(ids, ids)]
    q = len(ids)
    F = float(b @ np.linalg.inv(V) @ b) / q
    p = float(stats.f.sf(F, q, result["df"]))
    return {"F": F, "df_num": q, "df_den": result["df"], "p": p}


def linear_contrast(result, coefficients):
    c = np.zeros(len(result["names"]))
    for term, value in coefficients.items():
        c[result["names"].index(term)] = value
    est = float(c @ result["beta"])
    se = float(np.sqrt(c @ result["cov"] @ c))
    crit = stats.t.ppf(0.975, result["df"])
    t = est / se
    p = float(2 * stats.t.sf(abs(t), result["df"]))
    return {"beta": est, "se": se, "ci_low": est-crit*se, "ci_high": est+crit*se, "t": t, "p": p}
