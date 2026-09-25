"""Access to the long-period records downloaded by scripts/fetch_modes_data.py."""
import os
import numpy as np
from obspy import read, UTCDateTime
from obspy.geodetics import gps2dist_azimuth

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "seismograms")
TOHOKU = dict(time=UTCDateTime("2011-03-11T05:46:24"), lat=38.30, lon=142.37, name="Tohoku")
STATIONS = {"II_BFO": (48.331, 8.330), "IU_ANMO": (34.946, -106.457), "IU_CTAO": (-20.088, 146.254),
            "IU_KIP": (21.420, -158.011)}


def tohoku(station="II_BFO", components="Z", rotate=True):
    """Stream of acceleration records (m/s^2, 1 sample/s) starting one hour before the event.
    With rotate=True the horizontals are rotated to radial (R) and transverse (T)."""
    st = read(os.path.join(DATA, "tohoku_%s.mseed" % station))
    for tr in st:  # IU stations use LH1/LH2 for the (nominally N/E) horizontals
        tr.stats.channel = tr.stats.channel.replace("LH1", "LHN").replace("LH2", "LHE")
    if rotate and "T" in components or "R" in components:
        slat, slon = STATIONS[station]
        _, _, baz = gps2dist_azimuth(TOHOKU["lat"], TOHOKU["lon"], slat, slon)
        st.rotate("NE->RT", back_azimuth=baz)
    out = st.select(channel="LH[%s]" % components)
    for tr in out:
        tr.data = tr.data.astype(float)
    return out


def hours(tr):
    """Times in hours relative to the earthquake origin."""
    return (tr.times() - (TOHOKU["time"] - tr.stats.starttime)) / 3600.0


def amplitude_spectrum(tr, t_start_h=0.0, t_end_h=None):
    """Hann-windowed amplitude spectrum (frequency in mHz) of the record between the given hours."""
    h = hours(tr)
    sel = h >= t_start_h
    if t_end_h is not None:
        sel &= h <= t_end_h
    x = tr.data[sel] - tr.data[sel].mean()
    n = len(x)
    dt = tr.stats.delta
    f = np.fft.rfftfreq(n, dt) * 1e3
    return f, np.abs(np.fft.rfft(x * np.hanning(n))) * dt
