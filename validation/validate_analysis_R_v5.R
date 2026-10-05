#!/usr/bin/env Rscript

# Cross-software validation of the primary NHANES analysis
# ---------------------------------------------------------
# Purpose:
#   Independently reproduce the locked primary Python results using R's
#   established 'survey' package.
#
# IMPORTANT:
#   This is a validation implementation, not a translation of the Python
#   survey-estimation code. It starts from the raw NHANES XPT files and uses
#   survey::svydesign(), subset(), svyglm(), and regTermTest().
#
# Run:
#   Rscript validate_primary_analysis_R.R --data-dir /path/to/data_raw
#
# Required packages:
#   haven, survey
#
# Install once if needed:
#   install.packages(c("haven", "survey"))

options(stringsAsFactors = FALSE)
options(survey.lonely.psu = "adjust")

required <- c("haven", "survey")
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required R package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nInstall once with:\n  install.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(haven)
  library(survey)
})

# -----------------------------
# Command-line argument parsing
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || args[1] != "--data-dir") {
  stop(
    "Usage:\n  Rscript validate_primary_analysis_R.R --data-dir /path/to/data_raw",
    call. = FALSE
  )
}
data_dir <- normalizePath(args[2], mustWork = TRUE)

# -----------------------------
# Locked study specification
# -----------------------------
MIN_AGE <- 20
HBA1C_MAX <- 5.7
WEIGHT_TOL <- 1e-12
TOTAL_YEARS <- 5.2

expected <- list(
  n = 3079L,
  design_df = 40,
  interaction_F = 5.6863452144,
  interaction_p = 0.00670739084,
  slope_20_39 = 0.02662879194,
  slope_40_54 = 0.05112736866,
  slope_55_plus = -0.01890406625,
  ci20_low = 0.00380301065,
  ci20_high = 0.04945457323,
  p20 = 0.02336148342,
  ci40_low = 0.01324356846,
  ci40_high = 0.08901116886,
  p40 = 0.00942854864,
  ci55_low = -0.05132335076,
  ci55_high = 0.01351521826,
  p55 = 0.24555562443
)

# We expect close numerical agreement, not necessarily bit-for-bit identity.
# 1e-6 is stringent enough to detect a meaningful implementation discrepancy
# while allowing harmless floating-point/library differences.
TOL <- 1e-6

# -----------------------------
# Utilities
# -----------------------------
read_xpt_checked <- function(path) {
  if (!file.exists(path)) stop("Required NHANES file not found: ", path, call. = FALSE)
  d <- haven::read_xpt(path)
  if (!("SEQN" %in% names(d))) stop(basename(path), " is missing SEQN", call. = FALSE)
  if (anyDuplicated(d$SEQN)) stop(basename(path), " contains duplicate SEQN values", call. = FALSE)
  d
}

nhanes_filename <- function(key, prefix, suffix) {
  stems <- c(
    demo = "DEMO", bmx = "BMX", diq = "DIQ", fast = "FASTQX",
    ghb = "GHB", glu = "GLU", ins = "INS", hscrp = "HSCRP",
    tchol = "TCHOL", hdl = "HDL"
  )
  paste0(prefix, stems[[key]], suffix)
}

left_merge_one_to_one <- function(x, y, cols) {
  keep <- intersect(cols, names(y))
  y2 <- y[, keep, drop = FALSE]
  if (anyDuplicated(y2$SEQN)) stop("Duplicate SEQN in merge source", call. = FALSE)

  # Preserve DEMO row order exactly.
  idx <- match(x$SEQN, y2$SEQN)
  for (nm in setdiff(names(y2), "SEQN")) x[[nm]] <- y2[[nm]][idx]
  x
}

build_period <- function(label, prefix, suffix, fasting_weight, duration_years) {
  keys <- c("demo","bmx","diq","fast","ghb","glu","ins","hscrp","tchol","hdl")
  files <- setNames(
    lapply(keys, function(k) {
      read_xpt_checked(file.path(data_dir, nhanes_filename(k, prefix, suffix)))
    }),
    keys
  )

  demo_cols <- c("SEQN","RIDAGEYR","RIAGENDR","SDMVSTRA","SDMVPSU")
  d <- as.data.frame(files$demo[, demo_cols])

  selections <- list(
    bmx   = c("SEQN","BMXWAIST","BMXHT","BMXBMI","BMXHIP"),
    diq   = c("SEQN","DIQ010"),
    fast  = c("SEQN","PHAFSTHR","PHAFSTMN"),
    ghb   = c("SEQN","LBXGH"),
    glu   = c("SEQN", fasting_weight, "LBXGLU"),
    ins   = c("SEQN","LBXIN"),
    hscrp = c("SEQN","LBXHSCRP"),
    tchol = c("SEQN","LBXTC"),
    hdl   = c("SEQN","LBDHDD")
  )

  for (k in names(selections)) {
    d <- left_merge_one_to_one(d, as.data.frame(files[[k]]), selections[[k]])
  }

  d$period <- label
  d$raw_fasting_weight <- d[[fasting_weight]]
  d$analysis_weight <- d$raw_fasting_weight * (duration_years / TOTAL_YEARS)

  d$AGE <- d$RIDAGEYR
  d$is_male <- ifelse(d$RIAGENDR == 1, 1,
                      ifelse(d$RIAGENDR == 2, 0, NA_real_))
  d$waist <- d$BMXWAIST
  d$height <- d$BMXHT
  d$BMI <- d$BMXBMI

  if ("BMXHIP" %in% names(d)) {
    d$WHR <- d$BMXWAIST / d$BMXHIP
  } else {
    d$WHR <- NA_real_
  }

  d$waist_to_height <- d$waist / d$height
  d$hs_CRP <- d$LBXHSCRP
  d$log_hs_CRP <- ifelse(d$hs_CRP > 0, log(d$hs_CRP), NA_real_)
  d$HOMA_IR <- (d$LBXIN * d$LBXGLU) / 405
  d$Non_HDL <- d$LBXTC - d$LBDHDD
  d$fasting_hours <- d$PHAFSTHR + d$PHAFSTMN / 60

  # Unique masked design IDs across periods.
  d$stratum_u <- paste(d$period, d$SDMVSTRA, sep = ":")
  d$psu_u <- paste(d$period, d$SDMVPSU, sep = ":")
  d
}

cat("\n============================================================\n")
cat("R CROSS-SOFTWARE VALIDATION: PRIMARY NHANES ANALYSIS\n")
cat("============================================================\n\n")
cat("R version:      ", R.version.string, "\n", sep = "")
cat("survey version: ", as.character(packageVersion("survey")), "\n", sep = "")
cat("haven version:  ", as.character(packageVersion("haven")), "\n\n", sep = "")

# -----------------------------
# Build data from raw XPT files
# -----------------------------
d15 <- build_period(
  label = "2015-16",
  prefix = "",
  suffix = "_I.xpt",
  fasting_weight = "WTSAF2YR",
  duration_years = 2.0
)

d20 <- build_period(
  label = "2017-Mar2020",
  prefix = "P_",
  suffix = ".xpt",
  fasting_weight = "WTSAFPRP",
  duration_years = 3.2
)

# Harmonize period-specific columns before stacking.
# Some NHANES variables are available only in one period (notably BMXHIP in
# 2017-March 2020). These variables are not used in the primary model, but
# base R rbind() requires identical columns.
all_cols <- union(names(d15), names(d20))

add_missing_columns <- function(x, all_cols) {
  missing_cols <- setdiff(all_cols, names(x))
  for (nm in missing_cols) x[[nm]] <- NA
  x[, all_cols, drop = FALSE]
}

d15 <- add_missing_columns(d15, all_cols)
d20 <- add_missing_columns(d20, all_cols)

d <- rbind(d15, d20)

d$age_group <- cut(
  d$AGE,
  breaks = c(19, 39, 54, Inf),
  labels = c("20-39", "40-54", "55+"),
  right = TRUE
)

# Explicit reference levels are part of the locked model parameterization.
d$age_group <- factor(d$age_group, levels = c("20-39", "40-54", "55+"))
d$period <- factor(d$period, levels = c("2015-16", "2017-Mar2020"))

# -----------------------------
# Eligibility / participant flow
# -----------------------------
age_ok <- !is.na(d$AGE) & d$AGE >= MIN_AGE
diab_ok <- !is.na(d$DIQ010) & d$DIQ010 == 2
hba1c_ok <- !is.na(d$LBXGH) & d$LBXGH < HBA1C_MAX
fasting_w_ok <- !is.na(d$raw_fasting_weight) & d$raw_fasting_weight > WEIGHT_TOL
labs_ok <- !is.na(d$LBXGLU) & !is.na(d$LBXIN)
positive_ok <- !is.na(d$HOMA_IR) & !is.na(d$hs_CRP) & d$HOMA_IR > 0 & d$hs_CRP > 0
covariate_complete <- complete.cases(d[, c("LBXGH","waist","AGE","is_male","Non_HDL")])

m0 <- rep(TRUE, nrow(d))
m1 <- m0 & age_ok
m2 <- m1 & diab_ok
m3 <- m2 & hba1c_ok
m4 <- m3 & fasting_w_ok
m5 <- m4 & labs_ok
m6 <- m5 & positive_ok
m7 <- m6 & covariate_complete

flow <- data.frame(
  stage = c(
    "NHANES 2015-March 2020 participants",
    "Age >= 20",
    "DIQ010 == 2",
    "HbA1c < 5.7%",
    "Qualifying fasting-subsample weight",
    "Valid fasting glucose and insulin",
    "Positive HOMA-IR and hs-CRP",
    "Complete primary covariates"
  ),
  remaining_n = c(sum(m0),sum(m1),sum(m2),sum(m3),sum(m4),sum(m5),sum(m6),sum(m7))
)
flow$excluded_at_stage <- c(NA, -diff(flow$remaining_n))

cat("PARTICIPANT FLOW\n")
print(flow, row.names = FALSE)
cat("\n")

d$primary_domain <- m7

# -----------------------------
# Survey design
# -----------------------------
# NHANES fasting-subsample weights are defined only for participants selected
# into the fasting subsample. Therefore, construct the survey design on the
# qualifying fasting-subsample frame (finite, non-nominal-zero fasting weight),
# then estimate the primary analytic population as a domain within that design.
#
# This avoids inventing zero weights for people who were never members of the
# fasting subsample while retaining fasting-subsample participants who are
# outside the final analytic domain for correct domain variance estimation.
fasting_design_frame <- d[
  !is.na(d$raw_fasting_weight) &
    d$raw_fasting_weight > WEIGHT_TOL &
    !is.na(d$analysis_weight),
  ,
  drop = FALSE
]

fasting_design <- survey::svydesign(
  ids = ~psu_u,
  strata = ~stratum_u,
  weights = ~analysis_weight,
  data = fasting_design_frame,
  nest = TRUE
)

primary_design <- subset(fasting_design, primary_domain)

analytic_n <- sum(d$primary_domain)
design_df <- survey::degf(primary_design)

cat("DESIGN CHECKS\n")
cat("Analytic n: ", analytic_n, "\n", sep = "")
cat("Survey df:  ", design_df, "\n\n", sep = "")

# -----------------------------
# Locked primary model
# -----------------------------
# Python/patsy specification:
# log_hs_CRP ~ HOMA_IR * C(age_group)
#              + C(period) * (LBXGH + waist + AGE + is_male + Non_HDL)
#
# In R:
primary_formula <- log_hs_CRP ~
  HOMA_IR * age_group +
  period * (LBXGH + waist + AGE + is_male + Non_HDL)

fit <- survey::svyglm(
  primary_formula,
  design = primary_design,
  family = gaussian()
)

cat("PRIMARY MODEL COEFFICIENTS\n")
print(coef(summary(fit)))
cat("\n")

# -----------------------------
# Joint HOMA-IR x age-group test
# -----------------------------
interaction_test <- survey::regTermTest(
  fit,
  ~ HOMA_IR:age_group,
  method = "Wald"
)

interaction_F <- unname(interaction_test$Ftest)
interaction_p_regTermTest <- unname(interaction_test$p)

# regTermTest() uses a residual-denominator-df convention for this multi-
# parameter test. Preserve that output as a software-default diagnostic, but
# also evaluate the identical Wald F statistic using the prespecified NHANES
# survey design degrees of freedom used by the canonical analysis.
interaction_num_df <- 2
interaction_design_df <- survey::degf(primary_design)
interaction_p_design_df <- pf(
  interaction_F,
  df1 = interaction_num_df,
  df2 = interaction_design_df,
  lower.tail = FALSE
)

cat("JOINT HOMA-IR x AGE-GROUP TEST\n")
cat("R survey::regTermTest() default:\n")
print(interaction_test)
cat(sprintf(
  "\nSame Wald F evaluated with prespecified survey design df:\nF(%d,%d) = %.12f, p = %.12f\n\n",
  interaction_num_df, interaction_design_df,
  interaction_F, interaction_p_design_df
))

# -----------------------------
# Age-specific HOMA-IR slopes
# -----------------------------
b <- coef(fit)
V <- vcov(fit)

# Inspect names defensively rather than assuming hidden ordering.
needed_names <- c(
  "HOMA_IR",
  "HOMA_IR:age_group40-54",
  "HOMA_IR:age_group55+"
)
missing_names <- setdiff(needed_names, names(b))
if (length(missing_names) > 0) {
  stop(
    "Expected coefficient name(s) not found: ",
    paste(missing_names, collapse = ", "),
    "\nAvailable names:\n",
    paste(names(b), collapse = "\n"),
    call. = FALSE
  )
}

L_young <- setNames(rep(0, length(b)), names(b))
L_mid   <- L_young
L_old   <- L_young

L_young["HOMA_IR"] <- 1

L_mid["HOMA_IR"] <- 1
L_mid["HOMA_IR:age_group40-54"] <- 1

L_old["HOMA_IR"] <- 1
L_old["HOMA_IR:age_group55+"] <- 1

slope_stats <- function(L, label) {
  est <- sum(L * b)
  se <- sqrt(as.numeric(t(L) %*% V %*% L))
  df <- survey::degf(primary_design)
  tval <- est / se
  pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
  crit <- qt(0.975, df = df)

  data.frame(
    age_group = label,
    beta = est,
    SE = se,
    CI_low = est - crit * se,
    CI_high = est + crit * se,
    t = tval,
    df = df,
    p = pval,
    row.names = NULL
  )
}

slopes <- rbind(
  slope_stats(L_young, "20-39"),
  slope_stats(L_mid, "40-54"),
  slope_stats(L_old, "55+")
)

cat("AGE-SPECIFIC HOMA-IR SLOPES\n")
print(slopes, digits = 12, row.names = FALSE)
cat("\n")

# -----------------------------
# Validation against frozen Python results
# -----------------------------
check_num <- function(label, expected_value, obtained_value, tol = TOL) {
  diff <- abs(expected_value - obtained_value)
  pass <- is.finite(diff) && diff <= tol
  data.frame(
    quantity = label,
    expected = expected_value,
    R_obtained = obtained_value,
    abs_difference = diff,
    tolerance = tol,
    status = ifelse(pass, "PASS", "FAIL"),
    row.names = NULL
  )
}

checks <- rbind(
  check_num("Analytic n", expected$n, analytic_n, 0),
  check_num("Survey df", expected$design_df, design_df, 0),
  check_num("Interaction F", expected$interaction_F, interaction_F),
  check_num("Interaction p (design df)", expected$interaction_p, interaction_p_design_df),
  check_num("Slope 20-39", expected$slope_20_39, slopes$beta[slopes$age_group == "20-39"]),
  check_num("20-39 CI low", expected$ci20_low, slopes$CI_low[slopes$age_group == "20-39"]),
  check_num("20-39 CI high", expected$ci20_high, slopes$CI_high[slopes$age_group == "20-39"]),
  check_num("20-39 p", expected$p20, slopes$p[slopes$age_group == "20-39"]),
  check_num("Slope 40-54", expected$slope_40_54, slopes$beta[slopes$age_group == "40-54"]),
  check_num("40-54 CI low", expected$ci40_low, slopes$CI_low[slopes$age_group == "40-54"]),
  check_num("40-54 CI high", expected$ci40_high, slopes$CI_high[slopes$age_group == "40-54"]),
  check_num("40-54 p", expected$p40, slopes$p[slopes$age_group == "40-54"]),
  check_num("Slope 55+", expected$slope_55_plus, slopes$beta[slopes$age_group == "55+"]),
  check_num("55+ CI low", expected$ci55_low, slopes$CI_low[slopes$age_group == "55+"]),
  check_num("55+ CI high", expected$ci55_high, slopes$CI_high[slopes$age_group == "55+"]),
  check_num("55+ p", expected$p55, slopes$p[slopes$age_group == "55+"])
)

cat("SOFTWARE-DEFAULT DF DIAGNOSTIC\n")
cat(sprintf(
  "regTermTest default p = %.12f (its printed denominator df may differ)\n",
  interaction_p_regTermTest
))
cat(sprintf(
  "Design-df p          = %.12f using denominator df = %d\n\n",
  interaction_p_design_df, interaction_design_df
))

cat("PYTHON-vs-R LOCKED VALIDATION\n")
print(checks, digits = 12, row.names = FALSE)
cat("\n")

# ============================================================
# SECONDARY CROSS-SOFTWARE VALIDATION
# ============================================================

joint_test <- function(fit, design) {
  jt <- survey::regTermTest(fit, ~ HOMA_IR:age_group, method = "Wald")
  Fval <- unname(jt$Ftest)
  ddf <- survey::degf(design)
  list(
    test = jt,
    F = Fval,
    p_default = unname(jt$p),
    design_df = ddf,
    p_design = pf(Fval, 2, ddf, lower.tail = FALSE)
  )
}

age_slopes <- function(fit, design) {
  b <- coef(fit)
  V <- vcov(fit)
  ddf <- survey::degf(design)
  make_L <- function(extra = NULL) {
    L <- setNames(rep(0, length(b)), names(b))
    L["HOMA_IR"] <- 1
    if (!is.null(extra)) L[extra] <- 1
    L
  }
  calc <- function(L, label) {
    est <- sum(L * b)
    se <- sqrt(as.numeric(t(L) %*% V %*% L))
    tv <- est / se
    crit <- qt(.975, ddf)
    data.frame(
      age_group = label, beta = est, SE = se,
      CI_low = est - crit * se, CI_high = est + crit * se,
      t = tv, df = ddf,
      p = 2 * pt(abs(tv), ddf, lower.tail = FALSE)
    )
  }
  rbind(
    calc(make_L(), "20-39"),
    calc(make_L("HOMA_IR:age_group40-54"), "40-54"),
    calc(make_L("HOMA_IR:age_group55+"), "55+")
  )
}

cat("============================================================\n")
cat("SECONDARY VALIDATION 1: CRUDE MODEL\n")
cat("============================================================\n")

crude_fit <- survey::svyglm(
  log_hs_CRP ~ HOMA_IR * age_group + period,
  design = primary_design,
  family = gaussian()
)
crude_joint <- joint_test(crude_fit, primary_design)
crude_slopes <- age_slopes(crude_fit, primary_design)

print(crude_joint$test)
cat(sprintf(
  "\nSame Wald F with survey design df:\nF(2,%d) = %.12f, p = %.12f\n\n",
  crude_joint$design_df, crude_joint$F, crude_joint$p_design
))
cat("AGE-SPECIFIC CRUDE HOMA-IR SLOPES\n")
print(crude_slopes, digits = 12, row.names = FALSE)

crude_checks <- rbind(
  check_num("Crude interaction F", 1.70699985, crude_joint$F),
  check_num("Crude interaction p (design df)", 0.19435864, crude_joint$p_design),
  check_num("Crude slope 20-39", 0.13079796, crude_slopes$beta[1]),
  check_num("Crude slope 40-54", 0.15534031, crude_slopes$beta[2]),
  check_num("Crude slope 55+", 0.10172919, crude_slopes$beta[3])
)

cat("\nCRUDE MODEL LOCKED VALIDATION\n")
print(crude_checks, digits = 12, row.names = FALSE)
cat("\n")

cat("============================================================\n")
cat("SECONDARY VALIDATION 2: LOG-HOMA SENSITIVITY\n")
cat("============================================================\n")

# Preserve the same adjusted specification but replace raw HOMA-IR with
# log(HOMA-IR). The local variable is intentionally still named HOMA_IR so
# the joint-test and contrast helpers remain identical.
log_design <- primary_design
log_design$variables$HOMA_IR_raw <- log_design$variables$HOMA_IR
log_design$variables$HOMA_IR <- log(log_design$variables$HOMA_IR)

log_fit <- survey::svyglm(
  log_hs_CRP ~
    HOMA_IR * age_group +
    period * (LBXGH + waist + AGE + is_male + Non_HDL),
  design = log_design,
  family = gaussian()
)
log_joint <- joint_test(log_fit, log_design)
log_slopes <- age_slopes(log_fit, log_design)

print(log_joint$test)
cat(sprintf(
  "\nSame Wald F with survey design df:\nF(2,%d) = %.12f, p = %.12f\n\n",
  log_joint$design_df, log_joint$F, log_joint$p_design
))
cat("AGE-SPECIFIC log(HOMA-IR) SLOPES (diagnostic)\n")
print(log_slopes, digits = 12, row.names = FALSE)

# Only rounded locked values were retained for this sensitivity analysis,
# so use tolerances that reflect that recorded precision.
log_checks <- rbind(
  check_num("Log-HOMA interaction F", 1.86, log_joint$F, tol = 0.005),
  check_num("Log-HOMA interaction p (design df)", 0.169, log_joint$p_design, tol = 0.0005)
)

cat("\nLOG-HOMA SENSITIVITY LOCKED VALIDATION\n")
print(log_checks, digits = 12, row.names = FALSE)
cat("\n")

all_checks <- rbind(checks, crude_checks, log_checks)

if (all(all_checks$status == "PASS")) {
  cat("============================================================\n")
  cat("FULL CROSS-SOFTWARE VALIDATION: PASS\n")
  cat("Primary + crude + log-HOMA sensitivity reproduced.\n")
  cat("============================================================\n")
  quit(status = 0)
} else {
  cat("============================================================\n")
  cat("FULL CROSS-SOFTWARE VALIDATION: DISCREPANCY DETECTED\n")
  cat("Do not tune the R analysis to force agreement.\n")
  cat("Inspect the failed target(s) above before changing any code.\n")
  cat("============================================================\n")
  quit(status = 2)
}
