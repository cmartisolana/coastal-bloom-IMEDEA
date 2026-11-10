"""Reproduce all trend estimates and the three figures from the bundled data (no download needed).

Usage:  uv run python scripts/reproduce.py
Runs scripts/compute_trends.py, then scripts/plot_figures.py.
"""
from pathlib import Path
import os, subprocess, sys, tempfile

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'coastal-bloom-mpl'))
for script in ['compute_trends.py', 'plot_figures.py']:
    subprocess.run([sys.executable, str(ROOT / 'scripts' / script)], cwd=ROOT, check=True)
