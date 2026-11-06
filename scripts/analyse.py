"""Recalculate all plotted estimates, matching CIs and two-sided asymptotic P values.

Run from any directory. Scientific definitions preserve the submitted analysis.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import norm, theilslopes
from statsmodels.tsa.seasonal import STL

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from coastal_bloom.trends import monthly_trend
DATA = ROOT / 'data/processed'
OUT = ROOT / 'results/tables'
OUT.mkdir(parents=True, exist_ok=True)
REGIONS = ['SW_Africa', 'NA_east']
records, summaries, polarity, polarity_shares, ekman = [], [], [], [], []


def annual_trend(years, y, lags=2):
    """OLS with Bartlett HAC and n/(n-2) finite-sample scaling."""
    years, y = np.asarray(years, float), np.asarray(y, float)
    if len(y) != 18 or not np.isfinite(y).all():
        raise ValueError('Expected 18 complete annual observations')
    X = np.column_stack([np.ones(len(y)), years-years.mean()])
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    scores = X * (y-X@b)[:, None]
    meat = scores.T@scores
    for lag in range(1, lags+1):
        cross = scores[lag:].T@scores[:-lag]
        meat += (1-lag/(lags+1))*(cross+cross.T)
    inv = np.linalg.inv(X.T@X)
    cov = inv@meat@inv * len(y)/(len(y)-2)
    return b[1], np.sqrt(max(cov[1, 1], 0))


def infer(region, metric, dates, y, units, annual=False, factor=1, sen=np.nan):
    if annual:
        slope, se = annual_trend(dates, y)
    else:
        h = monthly_trend(dates, y, 12)
        slope, se = h['slope_per_year'], h['hac_standard_error']
    slope, se = slope*factor, se*factor
    p = float(2*norm.sf(abs(slope/se))) if se > 0 else (1. if slope == 0 else 0.)
    lo, hi = slope-1.96*se, slope+1.96*se
    records.append(dict(region=region, metric=metric, estimate=slope, ci95_low=lo,
                        ci95_high=hi, standard_error=se, p_value=p, units=units,
                        n=int(np.isfinite(y).sum()), hac_lags=2 if annual else 12,
                        model='OLS annual + HAC' if annual else 'OLS monthly + month indicators + HAC',
                        descriptive_sen=sen, significant_unadjusted=bool(lo>0 or hi<0)))
    return slope, lo, hi, p


def shares(values, area, years):
    # Equal monthly weighting is intentional: reproduces the proposal calculation.
    total = pd.DataFrame(values*area).groupby(years).sum()
    return 100*total.div(total.sum(axis=1), axis=0)


for tag in REGIONS:
    z = np.load(DATA / f'{tag}_ring_monthly.npz', allow_pickle=False)
    months = z['months'].astype(str)
    dates = months.astype('datetime64[M]')
    years = np.array([int(m[:4]) for m in months])
    area, edges = z['ring_area_km2'], z['ring_limit_km']
    if not np.array_equal(months, pd.period_range('2003-01', '2020-12', freq='M').astype(str)):
        raise ValueError('Unexpected month coverage')
    for pol in ['all', 'cyclonic', 'anticyclonic']:
        occ = z[f'occupancy_{pol}']
        if occ.shape != (216, 20) or not np.isfinite(occ).all() or np.any((occ<0)|(occ>1)):
            raise ValueError('Invalid occupancy')
        centroid = (occ*edges).sum(1)/occ.sum(1)
        stl = STL(centroid, period=12, seasonal=999, seasonal_deg=0, robust=True).fit().trend
        sen = theilslopes(stl, np.arange(216)/12)[0]
        b, lo, hi, p = infer(tag, f'eddy_{pol}|centroid', dates, centroid, 'km/year', sen=sen)
        np.savez_compressed(OUT/f'{tag}_{pol}_centroid_monthly.npz', months=months, centroid=centroid, stl=stl)
        pd.DataFrame({'month':months, 'centroid_km':centroid, 'stl_trend_km':stl}).to_csv(OUT/f'{tag}_{pol}_centroid_monthly.csv',index=False)
        sh = shares(occ, area, years)
        assert np.allclose(sh.sum(axis=1), 100)
        near = sh.iloc[:, :6].sum(axis=1)
        nb,nlo,nhi,npval = infer(tag,f'eddy_{pol}|share_0_60',sh.index,near,'percentage points/decade',True,10)
        row = dict(region=tag, polarity=pol, centroid_mean_km=centroid.mean(),sen_km_per_yr=sen,
                   ols_km_per_yr=b,lo=lo,hi=hi,p_value=p, near_share_pct=near.mean(),near_share_pp_dec=nb,near_lo=nlo,near_hi=nhi)
        polarity.append(row)
        per_ring=[]
        for k in range(20):
            label=f'{10*k}-{10*k+10}'
            rb,rl,rh,rp=infer(tag,f'eddy_{pol}|share_ring_{label}',sh.index,sh[k],'percentage points/decade',True,10)
            rr=dict(ring_km=label,slope_pp_per_decade=rb,lo=rl,hi=rh,p_value=rp,mean_share=sh[k].mean())
            per_ring.append(rr);polarity_shares.append(dict(region=tag,polarity=pol,**rr))
        sh.columns=[f'{10*k}-{10*k+10}' for k in range(20)]
        sh.assign(nearshore_0_60=near).to_csv(OUT/f'{tag}_{pol}_occupancy_share_annual.csv',index_label='year')
        pd.DataFrame(occ,columns=sh.columns).assign(month=months).to_csv(OUT/f'{tag}_{pol}_occupancy_monthly.csv',index=False)
        if pol=='all':
            pd.DataFrame(per_ring).to_csv(OUT/f'{tag}_occupancy_share_ring_trends.csv',index=False)
            pd.DataFrame(dict(month=months,centroid_km=centroid,stl_trend_km=stl)).to_csv(OUT/f'{tag}_occupancy_centroid_monthly.csv',index=False)
            summaries.append(dict(region=tag,centroid_mean_km=centroid.mean(),centroid_sen_km_per_yr=sen,
                                  centroid_ols_km_per_yr=b,centroid_lo=lo,centroid_hi=hi,centroid_p_value=p,
                                  nearshore_share_pp_per_decade=nb,nearshore_lo=nlo,nearshore_hi=nhi))
    z = np.load(DATA/f'{tag}_ekman_rings.npz',allow_pickle=False)
    assert np.array_equal(z['months'].astype(str),months)
    for key in ['ekman_xshelf_15m','total_xshelf_15m']:
        values=z[f'ring_{key}']*100
        pd.DataFrame(values,columns=[f'{10*k}-{10*k+10}' for k in range(20)]).assign(month=months).to_csv(OUT/f'{tag}_{key}_monthly_cm_s.csv',index=False)
        for k in range(20):
            if not np.isfinite(values[:,k]).all(): continue
            metric=f'{key}|ring_{10*k}-{10*k+10}'
            b,lo,hi,p=infer(tag,metric,dates,values[:,k],'cm/s/decade',factor=10)
            ekman.append(dict(region=tag,metric=metric,slope_per_decade=b,lo=lo,hi=hi,p_value=p))
        for label,(a,c) in {'0_60':(0,6),'60_200':(6,20),'30_60':(3,6)}.items():
            y=np.nansum(values[:,a:c]*area[a:c],axis=1)/area[a:c][np.isfinite(values[0,a:c])].sum()
            infer(tag,f'{key}|band_{label}',dates,y,'cm/s/decade',factor=10)
    div=z['ring_ekman_div_15m'][:,3:]*1e7
    keep=np.isfinite(div).all(axis=0)
    dist, ar=edges[3:][keep],area[3:][keep]
    sh=shares(np.clip(div[:,keep],0,None),ar,years)
    assert np.allclose(sh.sum(axis=1),100)
    near=sh.loc[:,dist<=60].sum(axis=1)
    for label,y in [('upwelling_share_30_60',near)]+[(f'share_ring_{int(d)-10}-{int(d)}',sh[i]) for i,d in enumerate(dist)]:
        metric=f'ekman_div_15m|{label}'
        b,lo,hi,p=infer(tag,metric,sh.index,y,'percentage points/decade',True,10)
        ekman.append(dict(region=tag,metric=metric,slope_per_decade=b,lo=lo,hi=hi,p_value=p))
    sh.columns=[f'{int(d)-10}-{int(d)}' for d in dist]
    sh.to_csv(OUT/f'{tag}_upwelling_share_annual.csv',index_label='year')
    pd.DataFrame(z['ring_ekman_div_15m'],columns=[f'{10*k}-{10*k+10}' for k in range(20)]).assign(month=months).to_csv(OUT/f'{tag}_divergence_monthly_s-1.csv',index=False)

pd.DataFrame(summaries).to_csv(OUT/'equivalent_metrics_summary.csv',index=False)
pd.DataFrame(polarity).to_csv(OUT/'polarity_summary.csv',index=False)
pd.DataFrame(polarity_shares).to_csv(OUT/'polarity_share_ring_trends.csv',index=False)
pd.DataFrame(ekman).to_csv(OUT/'ekman_trends.csv',index=False)
pd.DataFrame(records).to_csv(OUT/'all_trends.csv',index=False)
print(f'Recomputed {len(records)} trend estimates with matching CIs and P values.')
table = pd.DataFrame(records)
table[table.metric.str.contains('centroid|share_0_60|upwelling_share_30_60|band_30_60')].to_csv(OUT/'key_results.csv',index=False)
