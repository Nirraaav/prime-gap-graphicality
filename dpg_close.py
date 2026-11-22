"""DPG property across t in (30.001, 30.417), using Proposition 5.1.

Usage:  python3 dpg_close.py [t0 t1 step]      (default: 29.99 30.43 0.01)

Proposition 5.1: N(sigma, T) <= 1.67e11 T^{1.448(1-sigma)} for sigma in [0.994, 1] and
exp(6.7e12) < T <= exp(2e13). See bellotti_extend.py for the constant.
Prints the grid steps that cannot be certified (expected: none) and the margins at a few points.
"""
import sys
import dpg_improved as dp
import improved as im
from improved import mp

LOGT_MAX_EXT = mp.mpf("2e13")
SIG_EXT, C_EXT = mp.mpf("0.994"), mp.mpf("1.67e11")


def margin_ext(t, theta):
    """DPG margin with Proposition 5.1 near sigma = 1; None outside its range of T."""
    logx, logy, logT = dp.block_logT(t, theta)
    if not (im.LOGT_MAX_B < logT <= LOGT_MAX_EXT):
        return None
    zs = im.log_zero_sum(logy, logT, dp.c0, "B24", sig=SIG_EXT, cB=C_EXT)
    return dp.margin_from_zero_sum(logx, logT, zs)


if __name__ == "__main__":
    t0, t1, step = (map(mp.mpf, sys.argv[1:4]) if len(sys.argv) == 4
                    else (mp.mpf("29.99"), mp.mpf("30.43"), mp.mpf("0.01")))
    bad = dp.scan(t0, t1, step, extra={"EXT": margin_ext})
    print(f"[{mp.nstr(t0, 7)}, {mp.nstr(t1, 7)}] step {step}: uncertified steps {bad}")
    for t in ("30.0", "30.1", "30.2", "30.3", "30.416"):
        t = mp.mpf(t)
        ms = [margin_ext(t, th) for th in dp.thetas if dp.short_inputs(t, th)]
        ms = [m for m in ms if m is not None]
        print(f"t = {t}: best margin with Proposition 5.1 = {mp.nstr(min(ms), 6)}")
