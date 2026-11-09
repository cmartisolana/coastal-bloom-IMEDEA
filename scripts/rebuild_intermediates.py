"""Rebuild monthly eddy/current inputs and maps from bundled intermediate data.

Writes results/rebuilt only, then compares them with the bundled inputs in data/processed.
The global products do not need to be downloaded for this verification.
"""
from pathlib import Path
import os, subprocess, sys, tempfile
import numpy as np
root=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'coastal-bloom-mpl'))
def run(name,*args):
    subprocess.run([sys.executable,str(root/'processing'/name),*args],cwd=root,check=True)
checks=[]
for tag in ['SW_Africa','NA_east']:
    run('contour_geometry.py',tag)
    run('eddy_occupancy.py',tag)
    run('s19_occurrence_map.py',tag)
    run('s28_ekman_metrics.py',tag)
    for name,keys in [
        (f'{tag}_ring_monthly.npz',['ring_area_km2','ring_limit_km','occupancy_all','occupancy_cyclonic','occupancy_anticyclonic']),
        (f'occurrence_map_{tag}.npz',['lon','lat','freq_pct','dist_km']),
        (f'{tag}_ekman_rings.npz',['ring_ekman_xshelf_15m','ring_total_xshelf_15m','ring_ekman_div_15m','map_mean_ekman_xshelf_15m'])]:
        a=np.load(root/'results/rebuilt'/name,allow_pickle=True)
        b=np.load(root/'data/processed'/name)
        for k in keys: np.testing.assert_allclose(a[k],b[k],rtol=1e-7,atol=1e-10,equal_nan=True,err_msg=f'{name}:{k}')
        checks.append(f'{name}: all selected arrays match data/processed')
    print(f'{tag}: processing reproduced',flush=True)
(root/'results/rebuilt/check.txt').write_text('\n'.join(checks)+'\n')
