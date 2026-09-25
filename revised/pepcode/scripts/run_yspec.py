"""Generate the yspec input files used for the Lecture 21 synthetics (acceleration seismograms, six hours, no attenuation)
and run yspec without gravitation (0) and with full self-gravitation (2).  Usage:
    YSPEC=/path/to/yspec python scripts/run_yspec.py
The outputs are written to data/yspec/ (already present in the repository), and are read by
scripts/L9F2.py.  yspec: https://github.com/da380/YSpec (Al-Attar & Woodhouse 2008)."""
import os, subprocess, shutil, sys
here = os.path.dirname(os.path.abspath(__file__)); out = os.path.join(here, "..", "data", "yspec")
os.makedirs(out, exist_ok=True)
yspec = os.environ.get("YSPEC")
if not yspec:
    sys.exit("set YSPEC to the yspec binary")
model = os.path.join(os.path.dirname(yspec), "..", "..", "examples", "prem.200")

def write(g, name):
    params = [("prefix for output files", name), ("earth model", "prem.200"), ("attenuation switch", "0"), ("gravitation", str(g)),
              ("output", "2"), ("potential and tilt corrections", "0"), ("lmin", "0"), ("lmax", "750"), ("fmin (mHz)", "0.2"),
              ("fmax (mHz)", "50"), ("length of time series (min)", "360"), ("time step (sec)", "1.0"), ("f11", "0.3"),
              ("f12", "0.6"), ("f21", "45"), ("f22", "50"), ("source depth (km)", "20.0"), ("source latitude", "0.0"),
              ("source longitude", "0.0"), ("Mrr", "1.0e21"), ("Mrt", "0.0"), ("Mrp", "0.0"), ("Mtt", "-1.0e21"), ("Mtp", "0.0"),
              ("Mpp", "0.0"), ("receiver depth", "3"), ("number of receivers", "3")]
    lines = ["# yspec parameter file", "# generated for the PEP lecture notes"]
    for c, v in params:
        lines += ["", "# " + c, "  " + v]
    lines += ["", "# receiver latitudes and longitudes", "  0  30", "  0  90", "  0  150", ""]
    open(os.path.join(out, f"in_{name}.txt"), "w").write("\n".join(lines))

for g in [0, 2]:
    write(g, f"q{g}")
if not os.path.exists(os.path.join(out, "prem.200")):
    shutil.copy(model, out)
for g in [0, 2]:
    subprocess.run([yspec, f"in_q{g}.txt"], cwd=out, check=True)
