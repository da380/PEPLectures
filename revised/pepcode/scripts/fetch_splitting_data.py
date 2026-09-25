"""Twenty days of LHZ after the 2011 Tohoku earthquake at a few quiet stations, for the
splitting figure of Lecture 24 (0S2 / 0S3 singlets).  Stored as acceleration in
data/seismograms/tohoku20_<net>_<sta>.mseed."""
import os
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "seismograms")
c = Client("EARTHSCOPE"); t0 = UTCDateTime("2011-03-11T05:46:24")
for net, sta, loc in [("II", "BFO", "00"), ("IU", "ANMO", "00"), ("IU", "CTAO", "00"), ("II", "PFO", "00"),
                      ("IU", "MAJO", "00"), ("II", "ESK", "00"), ("IU", "TUC", "00"), ("II", "NNA", "00")]:
    try:
        st = c.get_waveforms(net, sta, loc, "LHZ", t0 - 3600, t0 + 20 * 86400, attach_response=True)
        st.merge(method=1, fill_value="interpolate")
        st.remove_response(output="ACC", pre_filt=(0.00005, 0.0001, 0.4, 0.45), taper=True, taper_fraction=0.02)
        for tr in st:
            tr.data = tr.data.astype("float32")
        st.write(os.path.join(out, f"tohoku20_{net}_{sta}.mseed"), format="MSEED")
        print(sta, [(tr.stats.channel, tr.stats.npts) for tr in st], flush=True)
    except Exception as e:
        print(sta, "failed:", e, flush=True)
