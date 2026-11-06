"""Check scientific outputs against archived calculations and data invariants."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1]
O=ROOT/'results/tables';REF=ROOT/'data/reference'
checks=[]
def compare(file, keys, columns):
    old=pd.read_csv(REF/file).set_index(keys)
    new=pd.read_csv(O/file).set_index(keys)
    common=new.index.intersection(old.index)
    assert len(common)==len(new),file
    for col in columns:
        np.testing.assert_allclose(new.loc[common,col],old.loc[common,col],rtol=1e-7,atol=1e-9,err_msg=f'{file}: {col}')
    checks.append(f'{file}: {len(common)} rows match archived calculations')
compare('equivalent_metrics_summary.csv','region',['centroid_mean_km','centroid_sen_km_per_yr','centroid_ols_km_per_yr','centroid_lo','centroid_hi','nearshore_share_pp_per_decade','nearshore_lo','nearshore_hi'])
compare('polarity_summary.csv',['region','polarity'],['sen_km_per_yr','ols_km_per_yr','lo','hi'])
compare('polarity_share_ring_trends.csv',['region','polarity','ring_km'],['slope_pp_per_decade','lo','hi'])
compare('ekman_trends.csv',['region','metric'],['slope_per_decade','lo','hi'])
for tag in ['SW_Africa','NA_east']:
    compare(f'{tag}_occupancy_share_ring_trends.csv','ring_km',['slope_pp_per_decade','lo','hi'])
t=pd.read_csv(O/'all_trends.csv')
assert not t.duplicated(['region','metric']).any()
assert t[['estimate','ci95_low','ci95_high','standard_error','p_value']].notna().all().all()
np.testing.assert_allclose(t.ci95_low,t.estimate-1.96*t.standard_error)
np.testing.assert_allclose(t.ci95_high,t.estimate+1.96*t.standard_error)
np.testing.assert_allclose(t.p_value,2*norm.sf(abs(t.estimate/t.standard_error)))
assert t.p_value.between(0,1).all()
for file in O.glob('*share_annual.csv'):
    d=pd.read_csv(file).set_index('year')
    np.testing.assert_allclose(d[[c for c in d if '-' in c]].sum(axis=1),100)
# Independent statsmodels covariance check on one monthly and one annual series.
import sys
sys.path.insert(0,str(ROOT/'src'))
from coastal_bloom.trends import monthly_trend
import statsmodels.api as sm
z=np.load(ROOT/'data/processed/NA_east_ring_monthly.npz')
m=z['months'].astype('datetime64[M]');y=(z['occupancy_all']*z['ring_limit_km']).sum(1)/z['occupancy_all'].sum(1)
mi=m.astype(int);X=np.column_stack([np.ones(len(y)),(mi-mi[0])/12,(mi[:,None]%12==np.arange(1,12)).astype(float)])
f=sm.OLS(y,X).fit(cov_type='HAC',cov_kwds={'maxlags':12,'use_correction':True})
h=monthly_trend(m,y,12)
np.testing.assert_allclose([h['slope_per_year'],h['hac_standard_error']],[f.params[1],f.bse[1]],rtol=1e-9)
a=pd.read_csv(O/'NA_east_all_occupancy_share_annual.csv');X=sm.add_constant(a.year-a.year.mean())
f=sm.OLS(a.nearshore_0_60,X).fit(cov_type='HAC',cov_kwds={'maxlags':2,'use_correction':True})
r=t[(t.region=='NA_east')&(t.metric=='eddy_all|share_0_60')].iloc[0]
np.testing.assert_allclose([r.estimate,r.standard_error],[10*f.params.iloc[1],10*f.bse.iloc[1]],rtol=1e-9)
checks.append(f'{len(t)} estimates: finite matching intervals and two-sided normal P values')
checks.append('All annual shares sum to 100%; monthly and annual HAC match statsmodels independently')
(ROOT/'results/validation.txt').write_text('\n'.join(checks)+'\n')
print('\n'.join(checks))
