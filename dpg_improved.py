"""DPG property (Theorem 1.3) with the block estimate of Proposition 4.2.

Usage:  python3 dpg_improved.py [c0 t0 t1 step]      (default: 53.989 26.2 30.0 0.02)

For each grid step [t, t2], look for theta in {0.626, 0.628, ..., 0.748} such that, at both ends,
  (1) the short-interval argument gives a prime in (x, x + x^theta]   (thresholds.margin_short_interval), and
  (2) S(delta) <= 0.2499 P for delta = P / (2 x^theta)                (dpg_margin),
using Lemma 2.6 (Bellotti 2024) while log T <= 6.7e12, and Lemma 2.7 (Bellotti 2023) or
Lemma 2.8 (Ramare) above. Steps that cannot be certified are printed.
Proposition 5.1 (the extension above log T = 6.7e12) is added in dpg_close.py.
"""
import sys
import improved as im
from improved import mp, LN, lse
import thresholds as th

c0 = im.C0_PUBLISHED
K = im.kernel_constant(a=mp.mpf("0.2"))[0]
B23C, B23B, B23D = mp.mpf("4.45e12"), mp.mpf("57.8875"), mp.mpf(10393) / 900
thetas = [mp.mpf(x) / 1000 for x in range(626, 750, 2)]


def log_zero_sum_B23(logy, logT, us=mp.mpf("0.0004")):
    """log of S_y / y^2 with Lemma 2.8 on [0.6, 1-us] and Lemma 2.7 on [1-us, 1-eta]."""
    eta = mp.e ** im.log_eta(logT, c0)
    lNT = logT + LN(logT)
    terms = [lNT - logy, LN(2 * logy * mp.mpf("0.1")) + lNT - mp.mpf("0.8") * logy]
    rate = 2 * logy - mp.mpf(8) / 3 * (LN(3) + logT) - 2 * LN(logT)
    assert rate > 0
    terms.append(LN(2 * logy * 965 / rate) + mp.mpf(8) / 3 * us * (LN(3) + logT)
                 + (3 + 2 * us) * LN(logT) - 2 * us * logy)
    terms.append(LN(mp.mpf("51.5")) + 2 * LN(logT) - 2 * us * logy)
    f = lambda u: -2 * u * logy + B23B * u ** mp.mpf(1.5) * logT
    r = 2 * logy - mp.mpf(1.5) * B23B * mp.sqrt(us) * logT
    assert r > 0
    terms.append(LN(2 * logy * B23C / r) + B23D * LN(logT) + f(eta))
    return lse(terms)


def block_logT(t, theta):
    """log T of Proposition 4.2 for the top block y = x/2 and h = P/(2 x^theta)."""
    logx = mp.e ** t
    logh = logx * (1 - theta) - LN(2) - mp.mpf("1e-9")
    logy = logx - LN(2)
    ly2 = logy + LN(mp.mpf("2.02"))
    logT = LN(mp.mpf("18.6") * 100 / mp.mpf("0.005")) + logy + LN(ly2) + LN(LN(ly2)) - logh
    return logx, logy, logT


def margin_from_zero_sum(logx, logT, zs):
    cL = 1 / ((1 - mp.mpf("2.01") / 100) * (mp.mpf("0.999") - mp.mpf("0.005")) ** 2)
    top = LN(cL * 4 * K * (logT + 1)) + zs
    mg = LN(logx / LN(2) + 1) - LN(25) - 2 * LN(logx)
    return lse([top, mg]) - LN(mp.mpf("0.2499"))


def dpg_margin(t, theta, dens):
    """log(bound for S(delta)/P) - log(0.2499); negative means the DPG inequality holds.
    Returns None if the density input `dens` does not apply at this height."""
    logx, logy, logT = block_logT(t, theta)
    if dens == "B24":
        if logT > im.LOGT_MAX_B:
            return None
        zs = im.log_zero_sum(logy, logT, c0, "B24")
    elif dens == "B23":
        if logT < im.LOGT_MAX_B:
            return None
        zs = log_zero_sum_B23(logy, logT)
    else:
        zs = im.log_zero_sum(logy, logT, c0, "RAM")
    return margin_from_zero_sum(logx, logT, zs)


def short_inputs(t, theta):
    """Density inputs for which the short-interval argument gives a prime in (x, x + x^theta]."""
    out = set()
    for d in ("B24", "B23", "RAM"):
        try:
            m, logT = th.margin_short_interval(t, theta=theta, c0=c0, dens=d)
        except AssertionError:
            continue
        if d == "B24" and logT > th.LOGT_MAX_B24:
            continue
        if d == "B23" and logT < th.LOGT_MAX_B24:
            continue
        if m < 0:
            out.add(d)
    return out


def dpg_inputs(t, theta, extra=()):
    """Density inputs (names) for which dpg_margin < 0 at (t, theta). `extra` maps name -> margin function."""
    out = set()
    funcs = {d: (lambda d=d: dpg_margin(t, theta, d)) for d in ("B24", "B23", "RAM")}
    funcs.update({name: (lambda f=f: f(t, theta)) for name, f in dict(extra).items()})
    for name, f in funcs.items():
        try:
            m = f()
        except AssertionError:
            continue
        if m is not None and m < 0:
            out.add(name)
    return out


def good(t, t2, extra=()):
    """Return theta such that, for one fixed density input in each step, both the DPG inequality and the
    short-interval inequality hold at t AND at t2.  All margins are nonincreasing in t for fixed theta and
    input (Lemma 7.1), and each input is valid on an interval of t, so they then hold on all of [t, t2]."""
    for thv in thetas:
        if (dpg_inputs(t, thv, extra) & dpg_inputs(t2, thv, extra)) and (short_inputs(t, thv) & short_inputs(t2, thv)):
            return thv
    return None


def scan(t0, t1, step, extra=()):
    """Steps [t, t+step] of [t0, t1] that are NOT certified (expected: none)."""
    bad, t = [], t0
    while t < t1 - step / 1000:
        t2 = min(t + step, t1)
        if good(t, t2, extra) is None:
            bad.append((mp.nstr(t, 7), mp.nstr(t2, 7)))
        t = t2
    return bad


if __name__ == "__main__":
    if len(sys.argv) == 5:
        c0 = mp.mpf(sys.argv[1])
        t0, t1, step = map(mp.mpf, sys.argv[2:5])
    else:
        t0, t1, step = mp.mpf("26.2"), mp.mpf("30.0"), mp.mpf("0.02")
    bad = scan(t0, t1, step)
    print(f"c0={c0}: [{mp.nstr(t0, 7)}, {mp.nstr(t1, 7)}] step {step}: uncertified steps {bad}")
