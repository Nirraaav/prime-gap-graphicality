"""DPG for t in (30.001, 30.417) using published inputs only, when the next gap is not too large.

For each t, find the largest theta such that S(P / (2 x^theta)) <= 0.2499 P follows from
Lemma 2.6 (Bellotti 2024, log T <= 6.7e12). Then Lemma 2.2 gives the DPG property for every n with
p_{n+1} - p_n <= x^theta. The margin gets worse as theta grows, so theta is found by bisection.
This is the partial statement at the end of the remark "Where Proposition 5.1 is needed".
"""
import dpg_improved as dp
from improved import mp


def ok(t, th):
    try:
        m = dp.dpg_margin(t, th, "B24")
    except AssertionError:
        return False
    return m is not None and m < 0


if __name__ == "__main__":
    for t in [mp.mpf(s) for s in ("30.001", "30.1", "30.2", "30.3", "30.417")]:
        lo, hi = mp.mpf("0.3"), mp.mpf("0.75")
        assert ok(t, lo)
        while hi - lo > mp.mpf("1e-6"):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if ok(t, mid) else (lo, mid)
        logx = mp.e ** t
        print(f"t={mp.nstr(t, 6)}  max theta={mp.nstr(lo, 7)}  log of largest next gap={mp.nstr(lo * logx, 8)}")
