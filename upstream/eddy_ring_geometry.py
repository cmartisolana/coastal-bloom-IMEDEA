"""Ocean area of each 10-km coastal ring inside a region's Table S1 box.

Ocean pixels: 0.02 deg grid, land from Natural Earth 1:10 m land polygons
(rasterised with matplotlib Agg, no antialiasing; lakes fall inside land
polygons and are therefore excluded, as in the manuscript's ocean polygon).
Distance: great-circle distance to the filtered NE coastline (k-d tree on the
sphere), as in eddy_buffer_trend_regions.py. Rings [0,10), ..., [190,200].
Writes data/processed/rings_<tag>.npz with ring areas (km^2) and the pixel
distance field (float16, NaN on land) for plotting.
"""
from pathlib import Path as FilePath
import os
os.chdir(FilePath(__file__).resolve().parents[1])
FilePath("results/rebuilt").mkdir(parents=True, exist_ok=True)
import argparse, json, sys
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.spatial import cKDTree
sys.path.insert(0, "scripts")
def xyz(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)], -1)
def chord_to_km(c):
    return 6371.0088*2*np.arcsin(np.clip(c/2,0,1))

ap = argparse.ArgumentParser()
ap.add_argument("--tag", required=True); ap.add_argument("--regions", default="config/regions.json")
ap.add_argument("--coast", default="data/processed/ne10m_coast_filtered_2km.npz"); ap.add_argument("--lines", default="data/processed/ne10m_land_coast_lines.npz")
ap.add_argument("--res", type=float, default=0.02); ap.add_argument("--max-km", type=float, default=200); ap.add_argument("--ring-km", type=float, default=10)
a = ap.parse_args()
B = json.load(open(a.regions))[a.tag]["box"]
nx, ny = int(round((B[1]-B[0])/a.res)), int(round((B[3]-B[2])/a.res))
lon = B[0] + a.res*(np.arange(nx)+.5); lat = B[2] + a.res*(np.arange(ny)+.5)
# rasterise land
L = np.load(a.lines)["land"]
rings = [r[~np.isnan(r[:, 0])] for r in np.split(L, np.flatnonzero(np.isnan(L[:, 0])))]
rings = [r for r in rings if len(r) > 2 and r[:, 0].max() >= B[0] and r[:, 0].min() <= B[1] and r[:, 1].max() >= B[2] and r[:, 1].min() <= B[3]]
fig = plt.figure(figsize=(nx/100, ny/100), dpi=100); ax = fig.add_axes([0, 0, 1, 1]); ax.set_axis_off()
ax.set_xlim(B[0], B[1]); ax.set_ylim(B[2], B[3])
ax.add_collection(PolyCollection(rings, facecolors="k", edgecolors="none", antialiaseds=False))
fig.canvas.draw(); img = np.asarray(fig.canvas.buffer_rgba())[::-1, :, 0]; plt.close(fig)
ny, nx = img.shape  # rendered size may differ by a pixel from the request; grid follows the image
dx, dy = (B[1]-B[0])/nx, (B[3]-B[2])/ny
lon = B[0] + dx*(np.arange(nx)+.5); lat = B[2] + dy*(np.arange(ny)+.5)
land = img < 128
# distance
c = np.load(a.coast)
padlat = a.max_km/111.0 + .5; padlon = a.max_km/(111.0*np.cos(np.radians(min(80, max(abs(B[2]), abs(B[3])))))) + .5
cs = (c["coast_lon"] >= B[0]-padlon) & (c["coast_lon"] <= B[1]+padlon) & (c["coast_lat"] >= B[2]-padlat) & (c["coast_lat"] <= B[3]+padlat)
tree = cKDTree(xyz(c["coast_lon"][cs], c["coast_lat"][cs]))
LON, LAT = np.meshgrid(lon, lat)
dist = np.full((ny, nx), np.nan)
oc = ~land
dist[oc] = chord_to_km(tree.query(xyz(LON[oc], LAT[oc]))[0])
edges = np.arange(0, a.max_km+a.ring_km/2, a.ring_km); nr = len(edges)-1
ring = np.clip(np.floor(dist/a.ring_km), 0, nr-1); ring[~(dist <= a.max_km)] = np.nan
pix_area = (dx*111.195)*(dy*111.195)*np.cos(np.radians(LAT))
area = np.array([pix_area[ring == k].sum() for k in range(nr)])
np.savez_compressed(f"data/processed/rings_{a.tag}.npz", box=np.array(B), res=a.res, lon=lon, lat=lat, edges=edges, area_km2=area,
                    dist=np.where(dist <= a.max_km, dist, np.nan).astype(np.float16),
                    note="ocean area per ring (Table S1 box, NE10m land rasterised, filtered NE coastline distance)")
print(a.tag, f"grid {nx}x{ny}, land frac {land.mean():.2f}, ring areas (10^3 km2):", np.round(area/1e3, 1))
