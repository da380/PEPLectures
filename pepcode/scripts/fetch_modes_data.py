"""Download long-period records of the 2011 Tohoku earthquake (Mw 9.1, 2011-03-11 05:46 UTC)
from the EarthScope FDSN service, remove the instrument response to acceleration, and store
them as MiniSEED in data/seismograms/.  Ten days of LHZ (and LHE/LHN) at a few quiet stations.
Usage:  poetry run python scripts/fetch_modes_data.py"""
import os
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "seismograms")
os.makedirs(out, exist_ok=True)
c = Client("EARTHSCOPE")
t0 = UTCDateTime("2011-03-11T05:46:24")
for net, sta, loc in [("II", "BFO", "00")]:   # only BFO is used (Lecture 23, Figs. 1-2)
    try:
        st = c.get_waveforms(net, sta, loc, "LH?", t0 - 3600, t0 + 10 * 86400, attach_response=True)
        st.merge(method=1, fill_value="interpolate")
        st.remove_response(output="ACC", pre_filt=(0.00005, 0.0001, 0.4, 0.45), taper=True, taper_fraction=0.02)
        for tr in st:
            tr.data = tr.data.astype("float32")
        st.write(os.path.join(out, f"tohoku_{net}_{sta}.mseed"), format="MSEED")
        print(sta, [(tr.stats.channel, tr.stats.npts) for tr in st])
    except Exception as e:
        print(sta, "failed:", e)
