"""Locations used by the figure scripts: figures are written to revised/figures."""
from pathlib import Path

FIG = Path(__file__).resolve().parents[2] / "figures"
FIG.mkdir(exist_ok=True)
(FIG / "anim").mkdir(exist_ok=True)
