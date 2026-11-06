"""Stage 2a: extract META4.0 observations for Benguela (SW_Africa) and Gulf Stream (NA_east).

Keeps every observation whose centre lies in the Table S1 box padded by PAD degrees, dated
2002-09-01..2020-12-31 (pre-history for the 2003 contact windows). Network semantics:
track = network id; trajectory key = (track, segment); track 0 = unlinked observations.
previous_obs/next_obs are converted to local indices (-1 if the linked obs is outside the
extraction or absent). Effective contours (30 points) are kept, decoded to degrees.
Output: interim/meta4_<tag>_<polarity>.npz
"""
from pathlib import Path as FilePath
import os
os.chdir(FilePath(__file__).resolve().parents[1])
FilePath("results/rebuilt").mkdir(parents=True, exist_ok=True)
import json, sys, time as _t, numpy as np, netCDF4 as nc
PAD = 3.0; T0 = np.datetime64("2002-09-01"); T1 = np.datetime64("2020-12-31")
regions = json.load(open("config/regions.json")); tags = ["SW_Africa", "NA_east"]
d0 = np.datetime64("1950-01-01")
FILES = ["data/raw/META4_DT_allsat_anticyclonic_19930101_20230908.nc", "data/raw/META4_DT_allsat_cyclonic_19930101_20230908.nc"]
for f in [x for x in FILES if not sys.argv[1:] or any(p in x for p in sys.argv[1:])]:
    st = _t.time(); pol = "anticyclonic" if "anticyclonic" in f else "cyclonic"
    d = nc.Dataset(f); d.set_auto_mask(False)
    lon = d["longitude"][:]; lon = np.where(lon > 180, lon-360, lon).astype("float32"); lat = d["latitude"][:]
    tm = d["time"][:].astype(np.int32)
    tsel = (tm >= (T0-d0).astype(int)) & (tm <= (T1-d0).astype(int))
    sel = {}
    for tag in tags:
        x0, x1, y0, y1 = regions[tag]["box"]
        sel[tag] = np.flatnonzero(tsel & (lon >= x0-PAD) & (lon <= x1+PAD) & (lat >= y0-PAD) & (lat <= y1+PAD))
        print(pol, tag, len(sel[tag]), "obs", flush=True)
    allidx = np.unique(np.concatenate(list(sel.values())))
    scal = {}
    for v in ["track", "segment", "previous_obs", "next_obs", "amplitude", "effective_radius", "speed_radius",
              "speed_average", "effective_area", "observation_number", "effective_contour_shape_error"]:
        scal[v] = d[v][:][allidx]
    # contours: chunked contiguous reads, keep selected rows
    clon = np.empty((len(allidx), 30), "float32"); clat = np.empty_like(clon)
    CH = 2_000_000; n = len(tm)
    for s in range(0, n, CH):
        e = min(n, s+CH); k = np.flatnonzero((allidx >= s) & (allidx < e))
        if not len(k): continue
        r = allidx[k]-s
        clon[k] = d["effective_contour_longitude"][s:e][r]; clat[k] = d["effective_contour_latitude"][s:e][r]
        print(f"  contours {e/n:.0%} {_t.time()-st:.0f}s", flush=True)
    clon = np.where(clon > 180, clon-360, clon)
    pos = {g: i for i, g in enumerate(allidx)}  # global -> position in allidx
    for tag in tags:
        idx = sel[tag]; loc = np.searchsorted(allidx, idx)
        g2l = -np.ones(len(allidx), np.int64); g2l[loc] = np.arange(len(idx))
        def to_local(gl):
            out = -np.ones(len(gl), np.int64); ok = gl >= 0
            p = np.searchsorted(allidx, gl[ok]); p = np.clip(p, 0, len(allidx)-1)
            hit = allidx[p] == gl[ok]; tmp = -np.ones(ok.sum(), np.int64); tmp[hit] = g2l[p[hit]]; out[ok] = tmp; return out
        rec = dict(time=tm[idx], lon=lon[idx], lat=lat[idx], global_index=idx,
                   prev_local=to_local(scal["previous_obs"][loc].astype(np.int64)), next_local=to_local(scal["next_obs"][loc].astype(np.int64)),
                   prev_global=scal["previous_obs"][loc], next_global=scal["next_obs"][loc],
                   contour_lon=clon[loc], contour_lat=clat[loc],
                   **{v: scal[v][loc] for v in scal if v not in ("previous_obs", "next_obs")})
        np.savez_compressed(f"data/intermediate/meta4_{tag}_{pol}.npz", pad=PAD, box=np.array(regions[tag]["box"]),
                            source=f, doi=d.doi, units=json.dumps({"time": "days since 1950-01-01", "amplitude": "m", "effective_radius": "m",
                            "speed_radius": "m", "speed_average": "m/s (circum-averaged rotational speed)", "effective_area": "m2"}), **rec)
        print(pol, tag, "saved", f"{_t.time()-st:.0f}s", "contour lon range", float(clon[loc].min()), float(clon[loc].max()), flush=True)
    d.close(); del lon, lat, tm, tsel, clon, clat, scal
