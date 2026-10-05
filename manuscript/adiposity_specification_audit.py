import tomllib, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0,'.')
from src.data import build_combined_dataset
from src.models import fit_with_covariates
from src.survey import survey_linear_regression, joint_wald_test, linear_contrast

cfg=tomllib.load(open('config/analysis_config.toml','rb'))
df, primary, flow=build_combined_dataset(Path('/mnt/data'),cfg)

# helper

def summarize(name, formula, res):
    terms=[n for n in res['names'] if n.startswith('HOMA_IR:C(age_group)')]
    test=joint_wald_test(res,terms)
    rows=[]
    for ag,extra in [('20-39',None),('40-54','HOMA_IR:C(age_group)[T.40-54]'),('55+','HOMA_IR:C(age_group)[T.55+]')]:
        c={'HOMA_IR':1.0}
        if extra: c[extra]=1.0
        z=linear_contrast(res,c)
        rows.append((ag,z['beta'],z['ci_low'],z['ci_high'],z['p']))
    out={'model':name,'n':res['n'],'df':res['df'],'F':test['F'],'p_interaction':test['p']}
    for ag,b,lo,hi,p in rows:
        tag={'20-39':'young','40-54':'mid','55+':'old'}[ag]
        out[f'beta_{tag}']=b; out[f'lo_{tag}']=lo; out[f'hi_{tag}']=hi; out[f'p_{tag}']=p
    return out

# 1 combined exact common sample for waist/BMI/WHtR
common = primary & df[['waist','BMI','waist_to_height']].notna().all(axis=1)
rows=[]
for adip in ['waist','BMI','waist_to_height']:
    cov=['LBXGH',adip,'AGE','is_male','Non_HDL']
    formula=f"log_hs_CRP ~ HOMA_IR * C(age_group) + C(period) * ({' + '.join(cov)})"
    res=survey_linear_regression(df,formula,common)
    rows.append(summarize(f'Combined common sample: {adip}',formula,res))

# Also each alternative on max available within primary (to match usual sensitivity style)
for adip in ['waist','BMI','waist_to_height']:
    dom=primary & df[adip].notna()
    cov=['LBXGH',adip,'AGE','is_male','Non_HDL']
    formula=f"log_hs_CRP ~ HOMA_IR * C(age_group) + C(period) * ({' + '.join(cov)})"
    res=survey_linear_regression(df,formula,dom)
    rows.append(summarize(f'Combined available sample: {adip}',formula,res))

# 2 P period common exact sample incl WHR. Keep full P survey design frame and domain scores zero outside analytic domain.
pmask=df['period'].eq('2017-Mar2020')
dp=df.loc[pmask].copy()
primP=primary.loc[pmask].copy()
commonP=primP & dp[['waist','BMI','waist_to_height','WHR']].notna().all(axis=1)
for adip in ['waist','BMI','waist_to_height','WHR']:
    cov=['LBXGH',adip,'AGE','is_male','Non_HDL']
    formula=f"log_hs_CRP ~ HOMA_IR * C(age_group) + {' + '.join(cov)}"
    res=survey_linear_regression(dp,formula,commonP)
    rows.append(summarize(f'2017-Mar2020 WHR-common sample: {adip}',formula,res))

# P available samples individually
for adip in ['waist','BMI','waist_to_height','WHR']:
    dom=primP & dp[adip].notna()
    cov=['LBXGH',adip,'AGE','is_male','Non_HDL']
    formula=f"log_hs_CRP ~ HOMA_IR * C(age_group) + {' + '.join(cov)}"
    res=survey_linear_regression(dp,formula,dom)
    rows.append(summarize(f'2017-Mar2020 available sample: {adip}',formula,res))

out=pd.DataFrame(rows)
out.to_csv('/mnt/data/adiposity_specification_audit.csv',index=False)
print(out.to_string(index=False))
print('\nCommon combined n', common.sum())
print('Common P n', commonP.sum())
print('P primary n', primP.sum())
