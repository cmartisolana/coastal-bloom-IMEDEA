"""Eddy occurrence frequency map: % of days 2003-2020 each 0.04-deg cell lies inside any META4.0 effective contour
(both polarities, union), on the 0-200 km Table S1 box domain (data/processed/rings_<tag>.npz). Usage: s19_occurrence_map.py TAG"""
from pathlib import Path as FilePath
import os
os.chdir(FilePath(__file__).resolve().parents[1])
FilePath("results/rebuilt").mkdir(parents=True, exist_ok=True)
import sys, numpy as np, pandas as pd
from matplotlib.path import Path
tag = sys.argv[1]; IM = "data/intermediate"
rz = np.load(f"data/processed/rings_{tag}.npz"); pdist = rz["dist"].astype("float64"); plon = rz["lon"]; plat = rz["lat"]
ny, nx = (pdist.shape[0]//2)*2, (pdist.shape[1]//2)*2
g4lon = plon[:nx].reshape(-1, 2).mean(1); g4lat = plat[:ny].reshape(-1, 2).mean(1); g4s = np.diff(g4lon).mean()
d4 = np.nanmean(pdist[:ny, :nx].reshape(ny//2, 2, nx//2, 2), axis=(1, 3))
d0 = np.datetime64("1950-01-01"); td = (np.arange(np.datetime64("2003-01-01"), np.datetime64("2021-01-01"))-d0).astype(int)
cnt = np.zeros((len(g4lat), len(g4lon))); P = {}
for pol in ["anticyclonic", "cyclonic"]:
    e = np.load(f"{IM}/meta4_{tag}_{pol}.npz"); t = e["time"].astype(int); L, A = e["contour_lon"], e["contour_lat"]
    touch = (L.max(1) > g4lon[0]) & (L.min(1) < g4lon[-1]) & (A.max(1) > g4lat[0]) & (A.min(1) < g4lat[-1]) & (t >= td[0]) & (t <= td[-1])
    P[pol] = (t, L, A, pd.Series(np.flatnonzero(touch)).groupby(t[touch]).apply(np.array).to_dict())
for tt in td:
    m = np.zeros(cnt.shape, bool)
    for pol, (t, L, A, byd) in P.items():
        for i in byd.get(tt, []):
            lo, la = L[i], A[i]
            i0 = max(0, int(np.floor((lo.min()-g4lon[0])/g4s))); i1 = min(len(g4lon), int(np.ceil((lo.max()-g4lon[0])/g4s))+1)
            j0 = max(0, int(np.floor((la.min()-g4lat[0])/g4s))); j1 = min(len(g4lat), int(np.ceil((la.max()-g4lat[0])/g4s))+1)
            if i1 <= i0 or j1 <= j0: continue
            GX, GY = np.meshgrid(g4lon[i0:i1], g4lat[j0:j1])
            m[j0:j1, i0:i1] |= Path(np.column_stack([lo, la])).contains_points(np.column_stack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
    cnt += m
freq = np.where(np.isfinite(d4), 100*cnt/len(td), np.nan)
np.savez_compressed(f"results/rebuilt/occurrence_map_{tag}.npz", lon=g4lon, lat=g4lat, freq_pct=freq, dist_km=d4, n_days=len(td))
print(tag, "domain-mean occurrence %.1f%%, max %.1f%%" % (np.nanmean(freq), np.nanmax(freq)))
