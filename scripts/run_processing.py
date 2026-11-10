"""Optional: recompute the figure inputs from the regional extracts of the source products.

Runs processing steps 1-4 for each region, writes to results/rebuilt/ and checks the results
against the bundled copies in data/processed/. No download is needed.

Usage:  uv run python scripts/run_processing.py
"""
from pathlib import Path
import os, subprocess, sys, tempfile
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'coastal-bloom-mpl'))
STEPS = ['step1_eddy_contour_distance.py', 'step2_eddy_occupancy.py',
         'step3_eddy_occurrence_map.py', 'step4_currents.py']
CHECKS = {
    'eddy_occupancy': ['ring_area_km2', 'ring_limit_km', 'occupancy_all', 'occupancy_cyclonic', 'occupancy_anticyclonic'],
    'eddy_occurrence_map': ['lon', 'lat', 'freq_pct', 'dist_km'],
    'currents': ['ring_ekman_xshelf_15m', 'ring_total_xshelf_15m', 'ring_ekman_div_15m', 'map_mean_ekman_xshelf_15m'],
}

report = []
for region in ['benguela', 'gulf_stream']:
    for step in STEPS:
        subprocess.run([sys.executable, str(ROOT / 'processing' / step), region], cwd=ROOT, check=True)
    for product, keys in CHECKS.items():
        name = f'{region}_{product}.npz'
        new = np.load(ROOT / 'results/rebuilt' / name, allow_pickle=True)
        bundled = np.load(ROOT / 'data/processed' / name)
        for k in keys:
            np.testing.assert_allclose(new[k], bundled[k], rtol=1e-7, atol=1e-10, equal_nan=True, err_msg=f'{name}:{k}')
        report.append(f'{name}: matches data/processed')
    print(f'{region}: processing reproduced', flush=True)
(ROOT / 'results/rebuilt/check.txt').write_text('\n'.join(report) + '\n')
