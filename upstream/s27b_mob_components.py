"""Monthly means of GlobCurrent velocity components (ue, ve at 15 m and 0 m; uo, vo at 15 m) for one region box,
2003-2020, read from the public ARCO zarr (same product/version as s27_mob_ekman.py). Projection on the offshore
normal and divergence are linear, so they are applied later to these monthly means on the device geometry."""
from pathlib import Path as FilePath
import os
os.chdir(FilePath(__file__).resolve().parents[1])
FilePath("results/rebuilt").mkdir(parents=True, exist_ok=True)
import sys, json, time, numpy as np, blosc, urllib.request, urllib.error
BASE = "https://s3.waw3-1.cloudferro.com/mdl-arco-geo-037/arco/MULTIOBS_GLO_PHY_MYNRT_015_003/cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m_202411/geoChunked.zarr"
tag = sys.argv[1]; box = {"NA_east": [-80.26, -60.86, 31.03, 46.09], "SW_Africa": [9.83, 19.77, -36.56, -13.27]}[tag]; out = sys.argv[2]
def get(path, tries=5):
    for k in range(tries):
        try:
            with urllib.request.urlopen(f"{BASE}/{path}", timeout=60) as r: return r.read()
        except urllib.error.HTTPError as e:
            if e.code in (403, 404): return None
            if k == tries-1: raise
        except Exception:
            if k == tries-1: raise
        time.sleep(2*(k+1))
meta = json.loads(get(".zmetadata"))["metadata"]
def coord(n): a = meta[f"{n}/.zarray"]; return np.frombuffer(blosc.decompress(get(f"{n}/0")), dtype=a["dtype"])[:a["shape"][0]].astype("float64")
t = coord("time"); lat = coord("latitude"); lon = coord("longitude"); elev = coord("elevation")
days = np.datetime64("1950-01-01T00") + (t*3600).astype("timedelta64[s]")
x0, x1, y0, y1 = box
ii = np.flatnonzero((lon >= x0-0.75) & (lon <= x1+0.75)); jj = np.flatnonzero((lat >= y0-0.75) & (lat <= y1+0.75))
i0, i1, j0, j1 = ii[0], ii[-1]+1, jj[0], jj[-1]+1
CT, CS = meta["ue/.zarray"]["chunks"][0], meta["ue/.zarray"]["chunks"][2]
tidx = np.flatnonzero((days >= np.datetime64("2003-01-01")) & (days < np.datetime64("2021-01-01")))
mon_all = days[tidx].astype("datetime64[M]"); um = np.unique(mon_all)
e15 = int(np.flatnonzero(elev == -15)[0]); e0 = int(np.flatnonzero(elev == 0)[0])
VARS = [("ue", e15, "ue_15m"), ("ve", e15, "ve_15m"), ("uo", e15, "uo_15m"), ("vo", e15, "vo_15m"), ("ue", e0, "ue_0m"), ("ve", e0, "ve_0m")]
sums = {lab: np.zeros((len(um), j1-j0, i1-i0)) for _, _, lab in VARS}; cnts = {lab: np.zeros((len(um), j1-j0, i1-i0)) for _, _, lab in VARS}
st = time.time()
for tc in sorted(set(tidx//CT)):
    kt = tidx[(tidx//CT) == tc]; ta, tb = kt[0]-tc*CT, kt[-1]-tc*CT+1; mi = np.searchsorted(um, days[kt].astype("datetime64[M]"))
    for var, e, lab in VARS:
        a = meta[f"{var}/.zarray"]; sf = meta[f"{var}/.zattrs"].get("scale_factor", 1.0)
        for cj in range(j0//CS, (j1-1)//CS+1):
            for ci in range(i0//CS, (i1-1)//CS+1):
                b = get(f"{var}/{tc}.{e}.{cj}.{ci}")
                if b is None: continue
                raw = np.frombuffer(blosc.decompress(b), dtype=a["dtype"]).reshape(CT, 1, CS, CS)[ta:tb, 0]
                v = np.where(raw == a["fill_value"], np.nan, raw*sf)
                sj, si = cj*CS, ci*CS; ys, ye = max(sj, j0), min(sj+CS, j1); xs, xe = max(si, i0), min(si+CS, i1)
                blk = v[:, ys-sj:ye-sj, xs-si:xe-si]; ok = np.isfinite(blk)
                np.add.at(sums[lab][:, ys-j0:ye-j0, xs-i0:xe-i0], mi, np.where(ok, blk, 0)); np.add.at(cnts[lab][:, ys-j0:ye-j0, xs-i0:xe-i0], mi, ok)
    print(tag, "time chunk", tc, f"{time.time()-st:.0f}s", flush=True)
res = {lab: (sums[lab]/np.where(cnts[lab] > 0, cnts[lab], np.nan)).astype("float32") for lab in sums}
np.savez_compressed(out, months=um.astype(str), lon=lon[i0:i1], lat=lat[j0:j1], product="MULTIOBS_GLO_PHY_MYNRT_015_003 cmems_obs-mob_glo_phy-cur_my_0.25deg_P1D-m_202411", **res)
print("saved", out, f"{time.time()-st:.0f}s")
