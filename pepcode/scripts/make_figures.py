"""Regenerate the Python-made figures: python scripts/make_figures.py [names...]
With no arguments every L*F*.py and S*F*.py script in this directory is run in order."""
import runpy, sys, glob, os, time
here = os.path.dirname(os.path.abspath(__file__))
names = sys.argv[1:] or sorted(os.path.basename(p)[:-3] for p in glob.glob(os.path.join(here, "[LS]*F*.py")))
for name in names:
    t0 = time.time()
    runpy.run_path(os.path.join(here, name + ".py"), run_name="__main__")
    print("%-8s done in %5.1f s" % (name, time.time() - t0))
