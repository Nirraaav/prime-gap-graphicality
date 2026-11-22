"""Locate the endpoints of the DPG ranges (Theorem 1.3) by bisection; certify_dpg.py then checks every step."""
import dpg_improved as dp
import dpg_close as dc
from improved import mp

EXT = {"EXT": dc.margin_ext}


def bisect(pred, a, b, want_first_true=True):
    """pred false at a, true at b (or the reverse if want_first_true=False); returns boundary to 1e-4"""
    while b - a > mp.mpf("1e-4"):
        m = (a + b) / 2
        ok = pred(m)
        if want_first_true:
            (a, b) = (m, b) if not ok else (a, m)
        else:
            (a, b) = (m, b) if ok else (a, m)
    return a, b


g_ext = lambda t: dp.good(t, t, EXT) is not None
g_pub = lambda t: dp.good(t, t) is not None
print("start (any inputs):", bisect(g_ext, mp.mpf("26.0"), mp.mpf("26.8")))
print("published only, end of first window:", bisect(g_pub, mp.mpf("29.9"), mp.mpf("30.2"), want_first_true=False))
print("published only, start of second window:", bisect(g_pub, mp.mpf("30.2"), mp.mpf("30.7")))
