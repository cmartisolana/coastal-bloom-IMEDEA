"""Processing step 2: monthly eddy occupancy in each 10 km coastal ring (all, cyclonic, anticyclonic eddies).

Usage:   python processing/step2_eddy_occupancy.py REGION
Inputs:  data/inputs/meta4_<region>_<polarity>.npz, data/processed/<region>_coastal_rings.npz,
         results/rebuilt/<region>_<polarity>_contour_distance.npz (step 1)
Output:  results/rebuilt/<region>_eddy_occupancy.npz (bundled copy: data/processed/)
"""
from pathlib import Path as FilePath
import os
os.chdir(FilePath(__file__).resolve().parents[1])
FilePath("results/rebuilt").mkdir(parents=True, exist_ok=True)
import sys, numpy as np, pandas as pd
from matplotlib.path import Path
R = 6371.0088; tag = sys.argv[1]; IM = "data/inputs"; OUT = "results/rebuilt"
os.makedirs(OUT, exist_ok=True)
edges = np.arange(0, 201, 10); dk = edges[1:].astype(float); nr = 20
M = pd.period_range("2003-01", "2020-12", freq="M").astype(str)
rz = np.load(f"data/processed/{tag}_coastal_rings.npz"); pdist = rz["dist"].astype("float64"); plon = rz["lon"]; plat = rz["lat"]; res = float(rz["res"])
ok = np.isfinite(pdist); ring = np.where(ok, np.minimum((np.nan_to_num(pdist)//10).astype(int), 19), -1)
pa = np.where(ok, (R*np.radians(res))**2*np.cos(np.radians(plat))[:, None], 0.0)
A_ring = np.bincount(ring[ok], weights=pa[ok], minlength=nr)
ring_m = {}
# Daily eddy occupancy by ring
d0 = np.datetime64("1950-01-01"); days = np.arange(np.datetime64("2003-01-01"), np.datetime64("2021-01-01"))
dmon = days.astype("datetime64[M]").astype(str); tt_days = (days-d0).astype(int)
ny, nx = (pdist.shape[0]//2)*2, (pdist.shape[1]//2)*2
sub_ring = ring[:ny, :nx].reshape(ny//2, 2, nx//2, 2).transpose(0, 2, 1, 3).reshape(ny//2, nx//2, 4)
sub_area = pa[:ny, :nx].reshape(ny//2, 2, nx//2, 2).transpose(0, 2, 1, 3).reshape(ny//2, nx//2, 4)
g4lon = plon[:nx].reshape(-1, 2).mean(1); g4lat = plat[:ny].reshape(-1, 2).mean(1); g4s = np.diff(g4lon).mean()
cats = ["all", "anticyclonic", "cyclonic"]
occ = {c: np.zeros((len(days), nr)) for c in cats}
P = {}
for pol in ["anticyclonic", "cyclonic"]:
    e = np.load(f"{IM}/meta4_{tag}_{pol}.npz"); g = np.load(f"results/rebuilt/{tag}_{pol}_contour_distance.npz")
    t = e["time"].astype(int)
    touch = np.isfinite(g["dmin_contour_inbox_km"]) & (g["dmin_contour_inbox_km"] <= 200)
    P[pol] = (t, e["contour_lon"], e["contour_lat"], touch)
by_day = {pol: pd.Series(np.arange(len(v[0]))).groupby(v[0]).apply(np.array).to_dict() for pol, v in P.items()}
for di, tt in enumerate(tt_days):
    masks = {c: np.zeros(sub_ring.shape[:2], bool) for c in cats}
    for pol, (t, L, A, touch) in P.items():
        for i in by_day[pol].get(tt, []):
            if not touch[i]: continue
            lo, la = L[i], A[i]
            i0 = max(0, int(np.floor((lo.min()-g4lon[0])/g4s))); i1 = min(len(g4lon), int(np.ceil((lo.max()-g4lon[0])/g4s))+1)
            j0 = max(0, int(np.floor((la.min()-g4lat[0])/g4s))); j1 = min(len(g4lat), int(np.ceil((la.max()-g4lat[0])/g4s))+1)
            if i1 <= i0 or j1 <= j0: continue
            GX, GY = np.meshgrid(g4lon[i0:i1], g4lat[j0:j1])
            ins = Path(np.column_stack([lo, la])).contains_points(np.column_stack([GX.ravel(), GY.ravel()])).reshape(GX.shape)
            for c in ("all", pol):
                masks[c][j0:j1, i0:i1] |= ins
    for c in cats:
        r_ = sub_ring[masks[c]].ravel(); a_ = sub_area[masks[c]].ravel(); kk = r_ >= 0
        occ[c][di] = np.bincount(r_[kk], weights=a_[kk], minlength=nr)/A_ring
mon_mean = lambda X: pd.DataFrame(X).groupby(dmon).mean().reindex(M).values
for c in cats: ring_m[f"occupancy_{c}"] = mon_mean(occ[c])
np.savez_compressed(f"{OUT}/{tag}_eddy_occupancy.npz", months=np.array(M), ring_area_km2=A_ring, ring_limit_km=dk, **ring_m)
