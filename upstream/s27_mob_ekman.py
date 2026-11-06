"""Read GlobCurrent (MULTIOBS_GLO_PHY_MYNRT_015_003, daily 0.25 deg, v202411) from the public ARCO zarr for one region
and one time chunk; compute per-cell daily cross-shelf velocities (Ekman 15 m and 0 m; total 15 m) and Ekman divergence
(15 m and 0 m), and save monthly cell means. Usage: s27_mob_ekman.py TAG coords | TAG CHUNK"""
from pathlib import Path as FilePath
import os
os.chdir(FilePath(__file__).resolve().parents[1])
FilePath("results/rebuilt").mkdir(parents=True, exist_ok=True)
import sys, os, json, time, numpy as np, blosc, urllib.request, urllib.error
from scipy.spatial import cKDTree
BASE = "https://s3.waw3-1.cloudferro.com/mdl-arco-geo-037/arco/MULTIOBS_GLO_PHY_MYNRT_015_003/cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m_202411/geoChunked.zarr"
IM = "data/intermediate"; R = 6371.0088
tag, mode = sys.argv[1], sys.argv[2]
def get(path, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(f"{BASE}/{path}", timeout=60) as r: return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404): return None      # chunk not stored (all fill, e.g. land)
            if k == tries-1: raise
        except Exception as e:
            if k == tries-1: raise
            time.sleep(2*(k+1))
meta = json.loads(get(".zmetadata"))["metadata"]
def coord(name):
    a = meta[f"{name}/.zarray"]; return np.frombuffer(blosc.decompress(get(f"{name}/0")), dtype=a["dtype"])[:a["shape"][0]]
cfile = f"{IM}/mob_coords.npz"
if not os.path.exists(cfile):
    np.savez(cfile, time=coord("time"), lat=coord("latitude"), lon=coord("longitude"), elev=coord("elevation"))
C = np.load(cfile); t_h = C["time"].astype("float64"); lat = C["lat"].astype("float64"); lon = C["lon"].astype("float64"); elev = C["elev"]
days = np.datetime64("1950-01-01T00") + (t_h*3600).astype("timedelta64[s]")
if mode == "coords":
    print("time", days[0], days[-1], len(days), "lat", lat[:2], "lon", lon[:2], "elev", elev); sys.exit()
box = json.load(open("config/regions.json"))[tag]["box"]; x0, x1, y0, y1 = box
ii = np.flatnonzero((lon >= x0-0.75) & (lon <= x1+0.75)); jj = np.flatnonzero((lat >= y0-0.75) & (lat <= y1+0.75))
i0, i1, j0, j1 = ii[0], ii[-1]+1, jj[0], jj[-1]+1
CH_T, CH_S = meta["ue/.zarray"]["chunks"][0], meta["ue/.zarray"]["chunks"][2]
tsel = (days >= np.datetime64("2003-01-01")) & (days < np.datetime64("2021-01-01")); tidx = np.flatnonzero(tsel)
tchunks = sorted(set(tidx//CH_T)); tc = tchunks[int(mode)]
kt = tidx[(tidx//CH_T) == tc]; ta, tb = kt[0]-tc*CH_T, kt[-1]-tc*CH_T+1
e15 = int(np.flatnonzero(elev == -15)[0]); e0 = int(np.flatnonzero(elev == 0)[0])
def read(var, e):
    out = np.full((tb-ta, j1-j0, i1-i0), np.nan, "float32"); a = meta[f"{var}/.zarray"]; sf = meta[f"{var}/.zattrs"].get("scale_factor", 1.0)
    for cj in range(j0//CH_S, (j1-1)//CH_S+1):
        for ci in range(i0//CH_S, (i1-1)//CH_S+1):
            b = get(f"{var}/{tc}.{e}.{cj}.{ci}")
            if b is None: continue
            raw = np.frombuffer(blosc.decompress(b), dtype=a["dtype"]).reshape(CH_T, 1, CH_S, CH_S)[ta:tb, 0]
            v = np.where(raw == a["fill_value"], np.nan, raw*sf).astype("float32")
            sj, si = cj*CH_S, ci*CH_S; ys, ye = max(sj, j0), min(sj+CH_S, j1); xs, xe = max(si, i0), min(si+CH_S, i1)
            out[:, ys-j0:ye-j0, xs-i0:xe-i0] = v[:, ys-sj:ye-sj, xs-si:xe-si]
    return out
st = time.time()
# offshore normals at grid points (smoothed coast, as Stage 1)
gl, ga = lon[i0:i1], lat[j0:j1]; LO, LA = np.meshgrid(gl, ga)
def xyz(lo, la):
    lo, la = np.radians(lo), np.radians(la); return np.stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)], -1)
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
ne_, nn_ = (v*E).sum(-1), (v*N).sum(-1); m_ = np.hypot(ne_, nn_); ne_ = (ne_/m_).reshape(LO.shape); nn_ = (nn_/m_).reshape(LO.shape)
dx = R*1e3*np.radians(0.25)*np.cos(np.radians(LA)); dy = R*1e3*np.radians(0.25)
def div(u, v):
    du = (np.roll(u, -1, 2)-np.roll(u, 1, 2))/(2*dx); dv = (np.roll(v, -1, 1)-np.roll(v, 1, 1))/(2*dy)
    d = du+dv; d[:, :, [0, -1]] = np.nan; d[:, [0, -1], :] = np.nan; return d     # NaN where a neighbour is land (propagated)
dd = days[tidx[(tidx//CH_T) == tc]]; mon = dd.astype("datetime64[M]"); um = np.unique(mon)
out = {}
for e, lab in [(e15, "15m"), (e0, "0m")]:
    ue, ve = read("ue", e), read("ve", e)
    out[f"ekman_xshelf_{lab}"] = ue*ne_+ve*nn_; out[f"ekman_div_{lab}"] = div(ue, ve)
    if lab == "15m":
        uo, vo = read("uo", e), read("vo", e); out["total_xshelf_15m"] = uo*ne_+vo*nn_; del uo, vo
    del ue, ve
res = {k: np.stack([np.nanmean(v[mon == m], 0) for m in um]).astype("float32") for k, v in out.items()}
np.savez_compressed(f"{IM}/mob_{tag}_chunk{mode}.npz", months=um.astype(str), lon=gl, lat=ga, n_east=ne_, n_north=nn_, **res)
print(tag, "chunk", mode, "of", len(tchunks), "days", len(dd), str(dd[0])[:10], str(dd[-1])[:10], f"{time.time()-st:.0f}s",
      "valid frac", round(float(np.isfinite(res["ekman_xshelf_15m"][0]).mean()), 2))
