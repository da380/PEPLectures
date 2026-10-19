"""Slides, Lecture 13 introduction: global seismicity, M >= 5 from 1976 to 2025 (USGS ComCat), coloured by
depth over NASA's Blue Marble (topography and bathymetry, December 2004), centred on the Pacific,
66 S to 84 N.
Catalogue and background are downloaded once and cached in data/seismicity/.  Written as PNG to
../figures/slides (the slides are made in LibreOffice).
Usage:  poetry run python scripts/slides_seismicity.py"""
import csv
import io
import os
import urllib.request
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.paths import FIG
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "seismicity")
os.makedirs(DATA, exist_ok=True)
CAT = os.path.join(DATA, "usgs_M5_1976-2025.csv")
BM = os.path.join(DATA, "world.topo.bathy.200412.3x5400x2700.jpg")
BM_URL = "https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73909/world.topo.bathy.200412.3x5400x2700.jpg"
LON0 = 150.0   # central longitude
LAT0, LAT1 = -66.0, 84.0   # latitude range shown (drops most of Antarctica)

if not os.path.exists(CAT):
    rows = []
    for year in range(1976, 2026):
        url = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&minmagnitude=5&orderby=time-asc"
               f"&starttime={year}-01-01&endtime={year + 1}-01-01")
        with urllib.request.urlopen(url) as r:
            for row in csv.DictReader(io.TextIOWrapper(r, encoding="utf-8")):
                rows.append((row["latitude"], row["longitude"], row["depth"], row["mag"]))
        print(year, len(rows))
    with open(CAT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["latitude", "longitude", "depth_km", "mag"])
        w.writerows(rows)
if not os.path.exists(BM):
    urllib.request.urlretrieve(BM_URL, BM)

lat, lon, dep, mag = np.loadtxt(CAT, delimiter=",", skiprows=1, unpack=True)
lon = (lon - LON0 + 180.0) % 360.0 - 180.0
img = plt.imread(BM)
nx = img.shape[1]
img = np.roll(img, -int(round((LON0 / 360.0) * nx)), axis=1)

order = np.argsort(dep)   # deep events on top
fig = plt.figure(figsize=(16, 16 * (LAT1 - LAT0) / 360.0), dpi=200)
ax = fig.add_axes([0, 0, 1, 1])
ax.imshow(0.45 * img / 255.0, extent=(-180, 180, -90, 90), interpolation="bilinear")
sc = ax.scatter(lon[order], lat[order], c=dep[order], s=0.6 * (mag[order] - 4.0) ** 2.5, cmap="plasma_r",
                vmin=0, vmax=650, lw=0, alpha=0.8)
ax.set_xlim(-180, 180); ax.set_ylim(LAT0, LAT1)
ax.set_axis_off()
cax = fig.add_axes([0.03, 0.08, 0.18, 0.022])
cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
cb.set_label("depth / km", color="w"); cb.ax.tick_params(colors="w"); cb.outline.set_edgecolor("w")
cb.solids.set_alpha(1)
(FIG / "slides").mkdir(exist_ok=True)
fig.savefig(str(FIG / "slides" / "seismicity_M5_1976-2025.png"))
print("events plotted:", lat.size)
