"""Cross-shelf Ekman and total (geostrophic + Ekman) metrics from GlobCurrent (config/stage9_polarity_ekman.json).
Inputs: interim/mob_SW_Africa_chunk{0..3}.npz (monthly cell means already projected, from s27) and/or
interim/mob_components_<tag>.npz (monthly component means, from mob_components.py; projected here).
Months split between two time chunks are combined weighted by their number of days.
Outputs: results/stage9/<tag>_ekman_*.csv|npz"""
from pathlib import Path as FilePath
import os
os.chdir(FilePath(__file__).resolve().parents[1])
FilePath("results/rebuilt").mkdir(parents=True, exist_ok=True)
import sys, glob, os, json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import theilslopes
from statsmodels.tsa.seasonal import STL
sys.path.insert(0, "src"); from coastal_bloom.trends import monthly_trend
IM = "data/intermediate"; O9 = "results/rebuilt"; os.makedirs(O9, exist_ok=True); R = 6371.0088
KEYS = ["ekman_xshelf_15m", "ekman_xshelf_0m", "total_xshelf_15m", "ekman_div_15m", "ekman_div_0m"]
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
    comp = f"{IM}/mob_components_{tag}.npz"
    if os.path.exists(comp):
        z = np.load(comp, allow_pickle=True); gl, ga = z["lon"], z["lat"]; box = json.load(open("config/regions.json"))[tag]["box"]
        ne_, nn_ = normals(gl, ga, box); LA = np.meshgrid(gl, ga)[1]
        F = {"ekman_xshelf_15m": z["ue_15m"]*ne_+z["ve_15m"]*nn_, "ekman_xshelf_0m": z["ue_0m"]*ne_+z["ve_0m"]*nn_,
             "total_xshelf_15m": z["uo_15m"]*ne_+z["vo_15m"]*nn_, "ekman_div_15m": div(z["ue_15m"], z["ve_15m"], LA), "ekman_div_0m": div(z["ue_0m"], z["ve_0m"], LA)}
        return z["months"].astype(str), gl, ga, F
    files = sorted(glob.glob(f"{IM}/mob_{tag}_chunk*.npz"), key=lambda f: int(f.split("chunk")[1].split(".")[0]))
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
def nw(x, y, lags=2):
    X = np.column_stack([np.ones(len(x)), x-x.mean()]); b = np.linalg.lstsq(X, y, rcond=None)[0]; e = y-X@b
    S = (X*e[:, None]).T@(X*e[:, None])
    for L in range(1, lags+1):
        G = (X[L:]*e[L:, None]).T@(X[:-L]*e[:-L, None]); S += (1-L/(lags+1))*(G+G.T)
    Bi = np.linalg.inv(X.T@X); V = Bi@S@Bi*len(y)/(len(y)-2); return b[1], np.sqrt(V[1, 1])
rows = []
for tag in sys.argv[1:] or ["SW_Africa", "NA_east"]:
    months, gl, ga, F = load(tag); mo = months.astype("datetime64[M]"); yr = np.array([int(m[:4]) for m in months])
    # ring weights: 0.02-deg ocean pixels (Table S1 box geometry) per 0.25-deg cell
    rz = np.load(f"data/processed/rings_{tag}.npz"); pdist = rz["dist"].astype("float64"); plon, plat = rz["lon"], rz["lat"]; res = float(rz["res"])
    ok = np.isfinite(pdist); ring = np.where(ok, np.minimum(np.nan_to_num(pdist)//10, 19), -1).astype(int)
    pa = np.where(ok, (R*np.radians(res))**2*np.cos(np.radians(plat))[:, None], 0.0)
    ci = np.rint((np.broadcast_to(plon, pdist.shape)-gl[0])/0.25).astype(int); cj = np.rint((np.broadcast_to(plat[:, None], pdist.shape)-ga[0])/0.25).astype(int)
    inb = ok & (ci >= 0) & (ci < len(gl)) & (cj >= 0) & (cj < len(ga))
    out = {"months": months}; prof = {}
    for k in KEYS:
        X = F[k]; good = np.isfinite(X).mean(0) > 0.95                              # time-consistent support
        u = inb & good[np.clip(cj, 0, len(ga)-1), np.clip(ci, 0, len(gl)-1)]
        W = np.zeros((len(ga)*len(gl), 20)); np.add.at(W, (cj[u]*len(gl)+ci[u], ring[u]), pa[u])
        Xf = np.nan_to_num(X.reshape(len(months), -1)); Rm = (Xf@W)/np.where(W.sum(0) > 0, W.sum(0), np.nan)
        out[f"ring_{k}"] = Rm; out[f"cov_{k}"] = W.sum(0)/np.bincount(ring[ok], weights=pa[ok], minlength=20)
        out[f"map_mean_{k}"] = np.nanmean(X, 0)
        prof[k] = Rm
    out["lon"], out["lat"] = gl, ga
    np.savez_compressed(f"{O9}/{tag}_ekman_rings.npz", **out)
    dk = np.arange(10, 201, 10.0); A = np.bincount(ring[ok], weights=pa[ok], minlength=20)
    for k in ["ekman_xshelf_15m", "ekman_xshelf_0m", "total_xshelf_15m"]:
        Rm = prof[k]*100                                                              # cm/s, offshore positive
        for lab, (r0, r1) in {"0_60": (0, 6), "60_200": (6, 20), "30_60": (3, 6)}.items():
            y = np.nansum(Rm[:, r0:r1]*A[r0:r1], 1)/A[r0:r1][np.isfinite(Rm[0, r0:r1])].sum(); h = monthly_trend(mo, y, 12)
            rows.append(dict(region=tag, metric=f"{k}|band_{lab}", mean=np.nanmean(y), slope_per_decade=10*h["slope_per_year"],
                             lo=10*(h["slope_per_year"]-1.96*h["hac_standard_error"]), hi=10*(h["slope_per_year"]+1.96*h["hac_standard_error"]), units="cm/s"))
        for r in range(20):
            if not np.isfinite(Rm[:, r]).all(): continue
            h = monthly_trend(mo, Rm[:, r], 12)
            rows.append(dict(region=tag, metric=f"{k}|ring_{10*r}-{10*r+10}", mean=Rm[:, r].mean(), slope_per_decade=10*h["slope_per_year"],
                             lo=10*(h["slope_per_year"]-1.96*h["hac_standard_error"]), hi=10*(h["slope_per_year"]+1.96*h["hac_standard_error"]), units="cm/s"))
    for k in ["ekman_div_15m", "ekman_div_0m"]:
        Dv = prof[k][:, 3:]*1e7; d3 = dk[3:]; A3 = A[3:]                           # 1e-7 s-1, rings 30-200 km
        keep = np.isfinite(Dv).all(0); Dv, d3, A3 = Dv[:, keep], d3[keep], A3[keep]
        wpos = np.clip(Dv, 0, None); C = (wpos*A3*d3).sum(1)/np.where((wpos*A3).sum(1) > 0, (wpos*A3).sum(1), np.nan)
        Ci = pd.Series(C).interpolate(limit_direction="both").values
        stl = STL(Ci, period=12, seasonal=999, seasonal_deg=0, robust=True).fit(); sen = theilslopes(stl.trend, np.arange(len(Ci))/12)[0]; h = monthly_trend(mo, C, 12)
        rows.append(dict(region=tag, metric=f"{k}|upwelling_centroid_30_200", mean=np.nanmean(C), slope_per_decade=10*sen, lo=10*(h["slope_per_year"]-1.96*h["hac_standard_error"]),
                         hi=10*(h["slope_per_year"]+1.96*h["hac_standard_error"]), units="km (Theil-Sen; CI from OLS)", note=f"OLS {10*h['slope_per_year']:.2f}; rings used {int(keep.sum())}; months without upwelling {int(np.isnan(C).sum())}"))
        out2 = {"centroid": C, "stl": stl.trend, "rings_km": d3}
        tot = pd.DataFrame(wpos*A3).groupby(yr).sum(); sh = 100*tot.div(tot.sum(1), axis=0)
        nb = d3 <= 60; near = sh.loc[:, nb].sum(1); b, se = nw(near.index.values.astype(float), near.values)
        rows.append(dict(region=tag, metric=f"{k}|upwelling_share_30_60", mean=near.mean(), slope_per_decade=10*b, lo=10*(b-1.96*se), hi=10*(b+1.96*se), units="pp"))
        for i, dd in enumerate(d3):
            b, se = nw(sh.index.values.astype(float), sh.iloc[:, i].values)
            rows.append(dict(region=tag, metric=f"{k}|share_ring_{int(dd)-10}-{int(dd)}", mean=sh.iloc[:, i].mean(), slope_per_decade=10*b, lo=10*(b-1.96*se), hi=10*(b+1.96*se), units="pp"))
        for i in range(Dv.shape[1]):
            h = monthly_trend(mo, Dv[:, i], 12)
            rows.append(dict(region=tag, metric=f"{k}|ring_{int(d3[i])-10}-{int(d3[i])}", mean=Dv[:, i].mean(), slope_per_decade=10*h["slope_per_year"],
                             lo=10*(h["slope_per_year"]-1.96*h["hac_standard_error"]), hi=10*(h["slope_per_year"]+1.96*h["hac_standard_error"]), units="1e-7 s-1"))
        np.savez_compressed(f"{O9}/{tag}_{k}_upwelling.npz", months=months, **out2)
T = pd.DataFrame(rows); fn = f"{O9}/ekman_trends.csv"
if os.path.exists(fn):
    old = pd.read_csv(fn); T = pd.concat([old[~old.region.isin(T.region.unique())], T], ignore_index=True)
T.to_csv(fn, index=False)
pd.set_option("display.width", 220)
print(T[T.metric.str.contains("band_|centroid|share_30")].round(3).to_string(index=False))
