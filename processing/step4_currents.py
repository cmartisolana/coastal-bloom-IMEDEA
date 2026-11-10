"""Processing step 4: offshore Ekman and total velocity and Ekman divergence at 15 m, averaged in each 10 km coastal ring.

Usage:   python processing/step4_currents.py REGION
Inputs:  data/inputs/globcurrent_<region>_*.npz (monthly GlobCurrent grids; Benguela in four time chunks whose
         split months are combined with day weights), data/processed/<region>_coastal_rings.npz,
         data/processed/ne10m_coast_filtered_2km.npz
Output:  results/rebuilt/<region>_currents.npz (bundled copy: data/processed/)
"""
from pathlib import Path as FilePath
import os
os.chdir(FilePath(__file__).resolve().parents[1])
FilePath("results/rebuilt").mkdir(parents=True, exist_ok=True)
import sys, glob, os, json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
IM = "data/inputs"; O9 = "results/rebuilt"; os.makedirs(O9, exist_ok=True); R = 6371.0088
KEYS = ["ekman_xshelf_15m", "total_xshelf_15m", "ekman_div_15m"]
CHUNK_DAYS = {0: ("2003-01-01", "2005-05-25"), 1: ("2005-05-26", "2011-08-06"), 2: ("2011-08-07", "2017-10-17"), 3: ("2017-10-18", "2020-12-31")}
def xyz(lo, la):
    lo, la = np.radians(lo), np.radians(la); return np.stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)], -1)
def normals(gl, ga, box):
    x0, x1, y0, y1 = box; LO, LA = np.meshgrid(gl, ga)
    cz = np.load("data/processed/ne10m_coast_filtered_2km.npz"); kl = ((cz["coast_lon"]+180) % 360)-180; ka = cz["coast_lat"]; grp = cz["group"]
    P = xyz(kl, ka); stp = R*np.arccos(np.clip((P[1:]*P[:-1]).sum(-1), -1, 1)); brk = np.flatnonzero((stp > 5) | (grp[1:] != grp[:-1]))+1
    Ps = np.empty_like(P); k13 = np.ones(13)/13
    for r in np.split(np.arange(len(kl)), brk):
        if len(r) < 13: Ps[r] = P[r].mean(0); continue
        for c in range(3): Ps[r, c] = np.convolve(np.pad(P[r, c], 6, mode="edge"), k13, "valid")
    Ps /= np.linalg.norm(Ps, axis=1, keepdims=True); near = (kl > x0-5) & (kl < x1+5) & (ka > y0-5) & (ka < y1+5)
    Q = xyz(LO.ravel(), LA.ravel()); j = cKDTree(Ps[near]).query(Q)[1]; v = Q-Ps[near][j]
    E = np.stack([-np.sin(np.radians(LO.ravel())), np.cos(np.radians(LO.ravel())), 0*LO.ravel()], -1)
    N = np.stack([-np.sin(np.radians(LA.ravel()))*np.cos(np.radians(LO.ravel())), -np.sin(np.radians(LA.ravel()))*np.sin(np.radians(LO.ravel())), np.cos(np.radians(LA.ravel()))], -1)
    ne_, nn_ = (v*E).sum(-1), (v*N).sum(-1); m_ = np.hypot(ne_, nn_); return (ne_/m_).reshape(LO.shape), (nn_/m_).reshape(LO.shape)
def div(u, v, LA):
    dx = R*1e3*np.radians(0.25)*np.cos(np.radians(LA)); dy = R*1e3*np.radians(0.25)
    d = (np.roll(u, -1, 2)-np.roll(u, 1, 2))/(2*dx) + (np.roll(v, -1, 1)-np.roll(v, 1, 1))/(2*dy)
    d[:, :, [0, -1]] = np.nan; d[:, [0, -1], :] = np.nan; return d
def load(tag):
    comp = f"{IM}/globcurrent_{tag}_monthly.npz"
    if os.path.exists(comp):
        z = np.load(comp, allow_pickle=True); gl, ga = z["lon"], z["lat"]; box = json.load(open("config/regions.json"))[tag]["box"]
        ne_, nn_ = normals(gl, ga, box); LA = np.meshgrid(gl, ga)[1]
        F = {"ekman_xshelf_15m": z["ue_15m"]*ne_+z["ve_15m"]*nn_,
             "total_xshelf_15m": z["uo_15m"]*ne_+z["vo_15m"]*nn_, "ekman_div_15m": div(z["ue_15m"], z["ve_15m"], LA)}
        return z["months"].astype(str), gl, ga, F
    files = sorted(glob.glob(f"{IM}/globcurrent_{tag}_chunk*.npz"), key=lambda f: int(f.split("chunk")[1].split(".")[0]))
    allm = pd.period_range("2003-01", "2020-12", freq="M").astype(str); acc = {k: 0.0 for k in KEYS}; wts = {k: 0.0 for k in KEYS}
    for f in files:
        c = int(f.split("chunk")[1].split(".")[0]); z = np.load(f, allow_pickle=True); gl, ga = z["lon"], z["lat"]
        d = pd.date_range(*CHUNK_DAYS[c], freq="D"); nd = pd.Series(1, index=d).groupby(d.to_period("M").astype(str)).sum()
        mi = [list(allm).index(m) for m in z["months"].astype(str)]; w = nd.reindex(z["months"].astype(str)).values.astype(float)
        for k in KEYS:
            if isinstance(acc[k], float): acc[k] = np.zeros((len(allm),)+z[k].shape[1:]); wts[k] = np.zeros_like(acc[k])
            x = z[k]; ok = np.isfinite(x); ww = w[:, None, None]*ok
            np.add.at(acc[k], mi, np.where(ok, x, 0)*w[:, None, None]); np.add.at(wts[k], mi, ww)
    F = {k: acc[k]/np.where(wts[k] > 0, wts[k], np.nan) for k in KEYS}
    return np.array(allm), gl, ga, F
for tag in sys.argv[1:] or ["benguela", "gulf_stream"]:
    months, gl, ga, F = load(tag)
    # ring weights: 0.02-deg ocean pixels (Table S1 box geometry) per 0.25-deg cell
    rz = np.load(f"data/processed/{tag}_coastal_rings.npz"); pdist = rz["dist"].astype("float64"); plon, plat = rz["lon"], rz["lat"]; res = float(rz["res"])
    ok = np.isfinite(pdist); ring = np.where(ok, np.minimum(np.nan_to_num(pdist)//10, 19), -1).astype(int)
    pa = np.where(ok, (R*np.radians(res))**2*np.cos(np.radians(plat))[:, None], 0.0)
    ci = np.rint((np.broadcast_to(plon, pdist.shape)-gl[0])/0.25).astype(int); cj = np.rint((np.broadcast_to(plat[:, None], pdist.shape)-ga[0])/0.25).astype(int)
    inb = ok & (ci >= 0) & (ci < len(gl)) & (cj >= 0) & (cj < len(ga))
    out = {"months": months}
    for k in KEYS:
        X = F[k]; good = np.isfinite(X).mean(0) > 0.95                              # time-consistent support
        u = inb & good[np.clip(cj, 0, len(ga)-1), np.clip(ci, 0, len(gl)-1)]
        W = np.zeros((len(ga)*len(gl), 20)); np.add.at(W, (cj[u]*len(gl)+ci[u], ring[u]), pa[u])
        Xf = np.nan_to_num(X.reshape(len(months), -1)); Rm = (Xf@W)/np.where(W.sum(0) > 0, W.sum(0), np.nan)
        out[f"ring_{k}"] = Rm; out[f"cov_{k}"] = W.sum(0)/np.bincount(ring[ok], weights=pa[ok], minlength=20)
        out[f"map_mean_{k}"] = np.nanmean(X, 0)
    out["lon"], out["lat"] = gl, ga
    np.savez_compressed(f"{O9}/{tag}_currents.npz", **out)
    print(tag, "rebuilt 15 m current and divergence ring fields")
