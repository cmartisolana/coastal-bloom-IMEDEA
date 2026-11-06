"""One-command reproduction from the bundled processed data."""
from pathlib import Path
import os, subprocess, sys, tempfile
root=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'coastal-bloom-mpl'))
for script in ['analyse.py','validate.py','plot_figures.py']:
    subprocess.run([sys.executable,str(root/'scripts'/script)],cwd=root,check=True)
