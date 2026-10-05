from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple
import numpy as np
import pandas as pd


def _read_xpt(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required NHANES file not found: {path}")
    df = pd.read_sas(path, format="xport")
    if "SEQN" not in df.columns:
        raise ValueError(f"{path.name} is missing SEQN")
    if df["SEQN"].duplicated().any():
        raise ValueError(f"{path.name} contains duplicate SEQN values")
    return df


def _filename(key: str, cfg: dict) -> str:
    stems = {
        "demo": "DEMO", "bmx": "BMX", "diq": "DIQ", "fast": "FASTQX",
        "ghb": "GHB", "glu": "GLU", "ins": "INS", "hscrp": "HSCRP",
        "tchol": "TCHOL", "hdl": "HDL",
    }
    return f"{cfg['prefix']}{stems[key]}{cfg['suffix']}"


def build_period(data_dir: Path, cfg: dict, total_years: float) -> pd.DataFrame:
    files: Dict[str, pd.DataFrame] = {
        k: _read_xpt(data_dir / _filename(k, cfg))
        for k in ["demo", "bmx", "diq", "fast", "ghb", "glu", "ins", "hscrp", "tchol", "hdl"]
    }

    # DEMO is the design frame. Analytic-domain score contributions will later be
    # set to zero outside the eligible fasting/complete-case domain.
    demo_cols = ["SEQN", "RIDAGEYR", "RIAGENDR", "SDMVSTRA", "SDMVPSU"]
    df = files["demo"][demo_cols].copy()

    selections = {
        "bmx": ["SEQN", "BMXWAIST", "BMXHT", "BMXBMI", "BMXHIP"],
        "diq": ["SEQN", "DIQ010"],
        "fast": ["SEQN", "PHAFSTHR", "PHAFSTMN"],
        "ghb": ["SEQN", "LBXGH"],
        "glu": ["SEQN", cfg["fasting_weight"], "LBXGLU"],
        "ins": ["SEQN", "LBXIN"],
        "hscrp": ["SEQN", "LBXHSCRP"],
        "tchol": ["SEQN", "LBXTC"],
        "hdl": ["SEQN", "LBDHDD"],
    }
    for key, wanted in selections.items():
        available = [c for c in wanted if c in files[key].columns]
        df = df.merge(files[key][available], on="SEQN", how="left", validate="one_to_one")

    df["period"] = cfg["label"]
    df["raw_fasting_weight"] = df[cfg["fasting_weight"]]
    df["analysis_weight"] = df["raw_fasting_weight"] * (float(cfg["duration_years"]) / total_years)

    df["AGE"] = df["RIDAGEYR"]
    df["is_male"] = df["RIAGENDR"].map({1.0: 1.0, 2.0: 0.0})
    df["waist"] = df.get("BMXWAIST")
    df["height"] = df.get("BMXHT")
    df["BMI"] = df.get("BMXBMI")
    df["WHR"] = np.nan
    if "BMXHIP" in df.columns:
        df["WHR"] = df["BMXWAIST"] / df["BMXHIP"]
    df["waist_to_height"] = df["waist"] / df["height"]
    df["hs_CRP"] = df["LBXHSCRP"]
    df["log_hs_CRP"] = np.where(df["hs_CRP"] > 0, np.log(df["hs_CRP"]), np.nan)
    df["HOMA_IR"] = (df["LBXIN"] * df["LBXGLU"]) / 405.0
    df["Non_HDL"] = df["LBXTC"] - df["LBDHDD"]
    df["fasting_hours"] = df["PHAFSTHR"] + df["PHAFSTMN"] / 60.0

    # Make public-use masked design identifiers unique across combined periods.
    df["stratum_u"] = df["period"].astype(str) + ":" + df["SDMVSTRA"].astype("Int64").astype(str)
    df["psu_u"] = df["period"].astype(str) + ":" + df["SDMVPSU"].astype("Int64").astype(str)
    return df


def build_combined_dataset(data_dir: Path, config: dict) -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    total_years = float(config["combined"]["total_years"])
    a = build_period(data_dir, config["period_2015_2016"], total_years)
    b = build_period(data_dir, config["period_2017_2020"], total_years)
    df = pd.concat([a, b], ignore_index=True, sort=False)

    young_max = int(config["age_groups"]["young_max"])
    mid_max = int(config["age_groups"]["midlife_max"])
    labels = config["age_groups"]["labels"]
    df["age_group"] = pd.cut(
        df["AGE"], bins=[19, young_max, mid_max, np.inf], labels=labels, right=True
    )

    s = config["study"]
    tol = float(s["fasting_weight_zero_tolerance"])
    age_ok = df["AGE"].ge(int(s["minimum_age"]))
    diab_ok = df["DIQ010"].eq(2)
    hba1c_ok = df["LBXGH"].lt(float(s["hba1c_upper_exclusive"]))
    fasting_w_ok = df["raw_fasting_weight"].fillna(0).gt(tol)
    labs_ok = df["LBXGLU"].notna() & df["LBXIN"].notna()
    positive_ok = df["HOMA_IR"].gt(0) & df["hs_CRP"].gt(0)
    covars = config["variables"]["primary_covariates"]
    covariate_complete = df[covars].notna().all(axis=1)

    primary = age_ok & diab_ok & hba1c_ok & fasting_w_ok & labs_ok & positive_ok & covariate_complete

    # Sequential participant-flow counts from the full combined NHANES participant frame.
    stages = []
    def stage(label: str, mask: pd.Series, previous: pd.Series | None = None):
        stages.append({
            "stage": label,
            "remaining_n": int(mask.sum()),
            "excluded_at_stage": None if previous is None else int(previous.sum() - mask.sum()),
        })

    m0 = pd.Series(True, index=df.index)
    stage("NHANES 2015-March 2020 participants", m0)
    m1 = m0 & age_ok; stage("Age >= 20", m1, m0)
    m2 = m1 & diab_ok; stage("DIQ010 == 2", m2, m1)
    m3 = m2 & hba1c_ok; stage("HbA1c < 5.7%", m3, m2)
    m4 = m3 & fasting_w_ok; stage("Qualifying fasting-subsample weight", m4, m3)
    m5 = m4 & labs_ok; stage("Valid fasting glucose and insulin", m5, m4)
    m6 = m5 & positive_ok; stage("Positive HOMA-IR and hs-CRP", m6, m5)
    m7 = m6 & covariate_complete; stage("Complete primary covariates", m7, m6)

    audit = pd.DataFrame(stages)
    return df, primary, audit
