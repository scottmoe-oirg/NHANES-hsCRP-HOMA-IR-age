#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import tomllib
import pandas as pd

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.data import build_combined_dataset
from src.models import fit_primary, fit_crude, fit_with_covariates


def load_config(path: Path):
    with open(path, "rb") as f:
        return tomllib.load(f)


def save_model(prefix: Path, formula, test, slopes):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    slopes.to_csv(prefix.with_name(prefix.name + "_slopes.csv"), index=False)
    pd.DataFrame([test]).to_csv(prefix.with_name(prefix.name + "_interaction.csv"), index=False)
    with open(prefix.with_name(prefix.name + "_formula.txt"), "w") as f:
        f.write(formula + "\n")


def validate_primary(config, domain, test, slopes):
    v = config["validation"]
    tol = float(v["absolute_tolerance"])
    observed = {
        "n": int(domain.sum()),
        "design_df": int(test["df_den"]),
        "interaction_F": float(test["F"]),
        "interaction_p": float(test["p"]),
        "slope_20_39": float(slopes.loc[slopes.age_group.eq("20-39"), "beta"].iloc[0]),
        "slope_40_54": float(slopes.loc[slopes.age_group.eq("40-54"), "beta"].iloc[0]),
        "slope_55_plus": float(slopes.loc[slopes.age_group.eq("55+"), "beta"].iloc[0]),
    }
    expected = {
        "n": int(v["expected_n"]), "design_df": int(v["expected_design_df"]),
        "interaction_F": float(v["expected_interaction_F"]),
        "interaction_p": float(v["expected_interaction_p"]),
        "slope_20_39": float(v["expected_slope_20_39"]),
        "slope_40_54": float(v["expected_slope_40_54"]),
        "slope_55_plus": float(v["expected_slope_55_plus"]),
    }
    checks = {}
    for k in observed:
        checks[k] = (observed[k] == expected[k]) if k in {"n", "design_df"} else abs(observed[k]-expected[k]) <= tol
    return observed, expected, checks


def audit_rows(model_name, test, slopes):
    row = {"model": model_name, "interaction_F": test["F"], "interaction_p": test["p"], "design_df": test["df_den"]}
    for ag in ["20-39", "40-54", "55+"]:
        s = slopes.loc[slopes.age_group.eq(ag)].iloc[0]
        tag = ag.replace("+", "plus").replace("-", "_")
        row[f"slope_{tag}"] = s.beta
        row[f"ci_low_{tag}"] = s.ci_low
        row[f"ci_high_{tag}"] = s.ci_high
        row[f"p_{tag}"] = s.p
    return row


def run(mode, data_dir: Path, out_dir: Path, config):
    df, domain, flow = build_combined_dataset(data_dir, config)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir/"qc").mkdir(exist_ok=True); (out_dir/"primary").mkdir(exist_ok=True); (out_dir/"audits").mkdir(exist_ok=True)

    flow.to_csv(out_dir/"qc"/"participant_flow.csv", index=False)
    analytic = df.loc[domain].copy()
    analytic[["SEQN","period","AGE","age_group","HOMA_IR","hs_CRP","log_hs_CRP","LBXGH","waist","height","BMI","WHR","waist_to_height","Non_HDL","is_male","fasting_hours","analysis_weight","stratum_u","psu_u"]].to_csv(out_dir/"qc"/"analytic_dataset_key_variables.csv", index=False)

    formula, res, test, slopes = fit_primary(df, domain, config)
    save_model(out_dir/"primary"/"primary", formula, test, slopes)
    observed, expected, checks = validate_primary(config, domain, test, slopes)
    validation = {"observed": observed, "expected": expected, "checks": checks, "all_passed": all(checks.values())}
    with open(out_dir/"qc"/"validation.json", "w") as f: json.dump(validation, f, indent=2)
    if not validation["all_passed"]:
        raise RuntimeError(f"Primary validation failed: {checks}")

    if mode in {"crude_adjusted", "full_audit"}:
        f0, r0, t0, s0 = fit_crude(df, domain)
        rows = [audit_rows("Crude", t0, s0), audit_rows("Fully adjusted", test, slopes)]
        pd.DataFrame(rows).to_csv(out_dir/"audits"/"crude_vs_adjusted.csv", index=False)
        save_model(out_dir/"audits"/"crude", f0, t0, s0)

    if mode in {"sequential_adjustment", "full_audit"}:
        covars = list(config["variables"]["primary_covariates"])
        rows = []
        _, _, t0, s0 = fit_crude(df, domain)
        rows.append(audit_rows("Crude", t0, s0))
        # Each covariate alone
        for c in covars:
            _, _, tt, ss = fit_with_covariates(df, domain, [c], period_specific=True)
            rows.append(audit_rows(f"+ {c} only", tt, ss))
        rows.append(audit_rows("Fully adjusted", test, slopes))
        pd.DataFrame(rows).to_csv(out_dir/"audits"/"each_covariate_alone.csv", index=False)

        # Drop one from full
        drop_rows = []
        for c in covars:
            keep = [x for x in covars if x != c]
            _, _, tt, ss = fit_with_covariates(df, domain, keep, period_specific=True)
            drop_rows.append(audit_rows(f"Full minus {c}", tt, ss))
        pd.DataFrame(drop_rows).to_csv(out_dir/"audits"/"drop_one_from_full.csv", index=False)

        # Pre-specified explanatory sequence starting with waist, then remaining covariates.
        sequence = ["waist", "LBXGH", "AGE", "is_male", "Non_HDL"]
        seq_rows = [audit_rows("Crude", t0, s0)]
        current = []
        for c in sequence:
            current.append(c)
            _, _, tt, ss = fit_with_covariates(df, domain, current, period_specific=True)
            seq_rows.append(audit_rows(" + ".join(current), tt, ss))
        pd.DataFrame(seq_rows).to_csv(out_dir/"audits"/"sequential_adjustment.csv", index=False)

    print("\nReproducible NHANES analysis completed")
    print("="*72)
    print(f"Analytic n: {domain.sum():,}")
    print(f"Design df: {test['df_den']}")
    print(f"Primary interaction: F({test['df_num']},{test['df_den']})={test['F']:.6f}, p={test['p']:.9f}")
    print(slopes[["age_group","beta","ci_low","ci_high","p"]].to_string(index=False))
    print(f"Validation: {'PASS' if validation['all_passed'] else 'FAIL'}")
    print(f"Outputs: {out_dir}")


def main():
    p = argparse.ArgumentParser(description="Reproduce the NHANES hsCRP-HOMA-IR primary analysis and audits.")
    p.add_argument("--config", default=str(ROOT/"config"/"analysis_config.toml"))
    p.add_argument("--data-dir", default=str(ROOT/"data_raw"), help="Directory containing NHANES XPT files")
    p.add_argument("--out-dir", default=str(ROOT/"outputs"))
    p.add_argument("--mode", choices=["primary","crude_adjusted","sequential_adjustment","full_audit"], default="full_audit")
    args = p.parse_args()
    run(args.mode, Path(args.data_dir), Path(args.out_dir), load_config(Path(args.config)))

if __name__ == "__main__":
    main()
