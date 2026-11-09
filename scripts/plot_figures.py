"""Reproduce the three figures (eddies, currents and Ekman, eddy polarity); inference labels use OLS with OLS CIs."""
import json, os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Polygon, Patch
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
IM = "data/processed"; O7 = O9 = "results/tables"
OUT = "figures"; os.makedirs(OUT, exist_ok=True)
INK, INK2, MUTED, GRID, NS, LAND = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#d3d2cd", "#dcdbd6"; S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"font.size": 6.5, "font.family": "DejaVu Sans", "axes.linewidth": 0.6, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.titlesize": 7, "axes.titleweight": "bold", "axes.titlelocation": "left",
                     "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "none"})
REG = [("SW_Africa", "Benguela"), ("NA_east", "Gulf Stream")]
boxes = json.load(open("config/regions.json"))
land = np.load("data/processed/ne10m_land_coast_lines.npz")["land"]; LANDR = np.split(land, np.flatnonzero(np.isnan(land[:, 0])))
def f2(v): return f"{0.0 if abs(v) < 0.005 else v:+.2f}"
def c2(v): return f"{v:+.3f}" if 0 < abs(v) < 0.005 else f"{v:+.2f}"
def f1(v): return f"{0.0 if abs(v) < 0.05 else v:+.1f}"
def sig(lo, hi): return (np.asarray(lo) > 0) | (np.asarray(hi) < 0)
def layout():
    fig = plt.figure(figsize=(7.09, 8.4))
    G = fig.add_gridspec(3, 1, height_ratios=[1.6, 1, 1], left=0.08, right=0.985, top=0.97, bottom=0.05, hspace=0.32)
    return fig, G[0].subgridspec(1, 2, width_ratios=[0.55, 1.25], wspace=0.15), G[1:].subgridspec(2, 2, wspace=0.28, hspace=0.42)
def basemap(ax, tag, lab, name):
    box = boxes[tag]["box"]
    for r in LANDR:
        r = r[np.isfinite(r[:, 0])]
        if len(r) > 2 and (np.ptp(r[:, 0]) > 0.3 or np.ptp(r[:, 1]) > 0.3) and r[:, 0].max() > box[0]-3 and r[:, 0].min() < box[1]+3 and r[:, 1].max() > box[2]-3 and r[:, 1].min() < box[3]+3:
            ax.add_patch(Polygon(r, closed=True, fc=LAND, ec=MUTED, lw=0.3, zorder=3))
    ax.set_xlim(box[0]-0.3, box[1]+0.3); ax.set_ylim(box[2]-0.3, box[3]+0.3); ax.set_aspect(1/np.cos(np.radians(np.mean(box[2:]))))
    ax.set_title(f"{lab}  {name}", pad=3); ax.locator_params(axis="x", nbins=4)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{abs(v):.0f}°{'E' if v >= 0 else 'W'}"))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{abs(v):.0f}°{'N' if v >= 0 else 'S'}"))
def colorbar(fig, ax, pm, label, extend="neither"):
    cax = ax.inset_axes([0.55, 0.14, 0.38, 0.035]); cb = fig.colorbar(pm, cax=cax, orientation="horizontal", extend=extend)
    cb.set_label(label, fontsize=6.5); cb.ax.tick_params(labelsize=6, width=0.5, length=2); cb.outline.set_linewidth(0.4)
def years(a):
    a.xaxis.set_major_locator(matplotlib.dates.YearLocator(4)); a.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
    a.grid(axis="y", color=GRID, lw=0.4)
def sharebars(a, x, s, lo, hi, w=8, off=0.0, color=None):
    s, lo, hi = map(np.asarray, (s, lo, hi)); c = np.where(sig(lo, hi), color if color else np.where(x < 60, S2, S1), NS)
    a.bar(x+off, s, width=w, color=c, lw=0)
    a.errorbar(x+off, s, yerr=[s-lo, hi-s], fmt="none", ecolor=INK2, elinewidth=0.5, capsize=0)
def ringaxis(a, ylabel, excluded=False, shade=False):
    a.axhline(0, color=INK, lw=0.5); a.axvline(60, color=INK2, lw=0.5, ls="--"); a.set_xlim(0, 200)
    if excluded or shade: a.axvspan(0, 30, color=GRID, lw=0, zorder=0)
    if excluded: a.text(15, 0.5, "excluded", transform=a.get_xaxis_transform(), rotation=90, ha="center", va="center", fontsize=6, color=INK2)
    a.set_xlabel("distance from coast (km)"); a.set_ylabel(ylabel); a.grid(axis="y", color=GRID, lw=0.4)
def headroom(a, f=0.22):
    lo_, hi_ = a.get_ylim(); a.set_ylim(lo_, hi_+f*(hi_-lo_))
def save(fig, name):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(f"{OUT}/{name}.{ext}", dpi=300)
    plt.close(fig)

# ======================= Eddies =======================
S = pd.read_csv(f"{O7}/equivalent_metrics_summary.csv").set_index("region"); P = pd.read_csv(f"{O9}/polarity_summary.csv")
fig, top, rows = layout()
for col, (tag, name) in enumerate(REG):
    s = S.loc[tag]; ax = fig.add_subplot(top[col])
    om = np.load(f"{IM}/occurrence_map_{tag}.npz"); gl, ga, F, D = om["lon"], om["lat"], om["freq_pct"], om["dist_km"]; st = np.diff(gl).mean()
    pm = ax.pcolormesh(gl-st/2, ga-st/2, F, vmin=0, vmax=40, cmap="Blues", shading="auto", rasterized=True)
    ax.contourf(gl, ga, np.where(np.isfinite(F), D, np.nan), levels=[0, 30], colors=["#ffffff"], alpha=0.55)
    ax.contour(gl, ga, D, levels=[60], colors=[INK2], linewidths=0.6, linestyles="--")
    ax.contour(gl, ga, D, levels=[s.centroid_mean_km], colors=[INK], linewidths=1.0)
    basemap(ax, tag, "ab"[col], name)
    if col == 1: colorbar(fig, ax, pm, "eddy occurrence (% of days)")
    a = fig.add_subplot(rows[0, col]); C = pd.read_csv(f"{O7}/{tag}_occupancy_centroid_monthly.csv"); t = pd.to_datetime(C.month); tt = (t-t.iloc[0]).dt.days/365.25
    a.plot(t, C.centroid_km, color=S1, lw=0.5, alpha=0.45)
    hs = [(a.plot(t, C.stl_trend_km, color=S1, lw=1.6)[0], f"all  {f2(s.centroid_ols_km_per_yr)} ({c2(s.centroid_lo)} to {c2(s.centroid_hi)})")]
    for pol, c, ls in [("cyclonic", INK, "-"), ("anticyclonic", S3, "-")]:
        z = np.load(f"{O9}/{tag}_{pol}_centroid_monthly.npz", allow_pickle=True); r = P[(P.region == tag) & (P.polarity == pol)].iloc[0]
        h, = a.plot(t, z["stl"], color=c, lw=0.9 if sig(r.lo, r.hi) else 0.7, ls=ls if sig(r.lo, r.hi) else (0, (3, 1.5)))
        hs.append((h, f"{pol}  {f2(r.ols_km_per_yr)} ({c2(r.lo)} to {c2(r.hi)})"))
    mid = C.stl_trend_km.mean()-s.centroid_sen_km_per_yr*tt.mean(); a.plot(t, mid+s.centroid_sen_km_per_yr*tt, color=S1, lw=0.7, ls=(0, (4, 2)))
    years(a); a.set_ylabel("eddy centroid (km)"); headroom(a, 0.42)
    leg = a.legend(*zip(*hs), loc="upper left", fontsize=5.6, handlelength=1.6, title="OLS km yr⁻¹ (95% CI); curves: STL", title_fontsize=5.6)
    leg._legend_box.align = "left"; leg.get_title().set_color(INK2)
    a.set_title(f"{'cd'[col]}  {name}: cross-shore centroid", pad=3)
    a = fig.add_subplot(rows[1, col]); R = pd.read_csv(f"{O7}/{tag}_occupancy_share_ring_trends.csv"); x = np.arange(5, 200, 10)
    sharebars(a, x, R.slope_pp_per_decade, R.lo, R.hi); ringaxis(a, "share trend (pp decade⁻¹)", shade=True); headroom(a)
    a.text(0.02, 0.97, f"0–60 km {f1(s.nearshore_share_pp_per_decade)} pp decade⁻¹ ({f1(s.nearshore_lo)} to {f1(s.nearshore_hi)})", transform=a.transAxes, ha="left", va="top", fontsize=6, color=INK2)
    a.set_title(f"{'ef'[col]}  {name}: contribution share by ring", pad=3)
save(fig, "eddies")

# ======================= Currents and Ekman =======================
K = pd.read_csv(f"{O9}/ekman_trends.csv")
def krows(tag, prefix):
    s = K[(K.region == tag) & K.metric.str.startswith(prefix)].copy(); s["x"] = [int(m.split("_")[-1].split("-")[0])+5 for m in s.metric]; return s.sort_values("x")
fig, top, rows = layout()
for col, (tag, name) in enumerate(REG):
    z = np.load(f"{IM}/{tag}_ekman_rings.npz", allow_pickle=True); gl, ga = z["lon"], z["lat"]; rz = np.load(f"data/processed/rings_{tag}.npz"); pdist = rz["dist"].astype(float)
    ax = fig.add_subplot(top[col])
    LO, LA = np.meshgrid(gl, ga); ii = np.clip(np.rint((LO-rz["lon"][0])/float(rz["res"])).astype(int), 0, len(rz["lon"])-1)
    jj = np.clip(np.rint((LA-rz["lat"][0])/float(rz["res"])).astype(int), 0, len(rz["lat"])-1); inside = np.isfinite(pdist[jj, ii])
    pm = ax.pcolormesh(gl-0.125, ga-0.125, np.where(inside, 100*z["map_mean_ekman_xshelf_15m"], np.nan), norm=TwoSlopeNorm(0, -6, 6), cmap="RdBu_r", shading="auto", rasterized=True, lw=0)
    ax.contour(rz["lon"], rz["lat"], pdist, levels=[60], colors=[INK2], linewidths=0.6, linestyles="--")
    basemap(ax, tag, "ab"[col], name)
    if col == 1: colorbar(fig, ax, pm, "offshore Ekman velocity (cm s⁻¹)", "both")
    a = fig.add_subplot(rows[0, col])
    for i, (key, c) in enumerate([("ekman_xshelf_15m|ring_", S3), ("total_xshelf_15m|ring_", INK2)]):
        s = krows(tag, key); sharebars(a, s.x.values.astype(float), s.slope_per_decade, s.lo, s.hi, w=4, off=-2+4*i, color=c)
    ringaxis(a, "trend (cm s⁻¹ decade⁻¹)", shade=True); headroom(a, 0.3)
    a.legend(handles=[Patch(color=S3, label="Ekman"), Patch(color=INK2, label="geostrophic + Ekman")], loc="upper right", fontsize=5.8, ncol=2, handlelength=1.2)
    a.set_title(f"{'cd'[col]}  {name}: offshore velocity at 15 m", pad=3)
    a = fig.add_subplot(rows[1, col]); s = krows(tag, "ekman_div_15m|share_ring_")
    sharebars(a, s.x.values.astype(float), s.slope_per_decade, s.lo, s.hi); ringaxis(a, "share trend (pp decade⁻¹)", excluded=True); headroom(a)
    b0 = K[(K.region == tag) & (K.metric == "ekman_div_15m|upwelling_share_30_60")].iloc[0]
    a.text(0.98, 0.97, f"30–60 km {f1(b0.slope_per_decade)} pp decade⁻¹ ({f1(b0.lo)} to {f1(b0.hi)})", transform=a.transAxes, ha="right", va="top", fontsize=6, color=INK2)
    a.set_title(f"{'ef'[col]}  {name}: positive-divergence share", pad=3)
save(fig, "currents_Ekman")

# ======================= Eddy polarity share =======================
PS = pd.read_csv(f"{O9}/polarity_share_ring_trends.csv")
fig, ax = plt.subplots(1, 2, figsize=(7.09, 2.4), gridspec_kw=dict(left=0.08, right=0.985, top=0.9, bottom=0.17, wspace=0.28))
for col, (tag, name) in enumerate(REG):
    a = ax[col]; x = np.arange(5, 200, 10)
    for i, (pol, c) in enumerate([("cyclonic", INK), ("anticyclonic", S3)]):
        s = PS[(PS.region == tag) & (PS.polarity == pol)]; sharebars(a, x, s.slope_pp_per_decade, s.lo, s.hi, w=4, off=-2+4*i, color=c)
    ringaxis(a, "share trend (pp decade⁻¹)", shade=True); headroom(a, 0.25)
    a.legend(handles=[Patch(color=INK, label="cyclonic"), Patch(color=S3, label="anticyclonic")], loc="upper left", fontsize=5.8, ncol=2, handlelength=1.2)
    a.set_title(f"{'ab'[col]}  {name}: contribution share by polarity", pad=3)
save(fig, "eddy_polarity_share"); print("saved", sorted(os.listdir(OUT)))
