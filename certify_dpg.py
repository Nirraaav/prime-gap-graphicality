"""Certification of Theorem 1.3 (DPG property for log log p_n >= 26.246).

Every grid step [t, t2] is certified with one fixed (theta, density input) valid at both ends
(see dpg_improved.good). Beyond t = 80 the Lemma 2.7 / Lemma 2.8 inputs remain valid and all margins
keep decreasing (Lemma 7.1), so the last certified step extends to infinity.

  with Proposition 5.1:    [26.246, 80]
  published inputs only:   [26.246, 30.001] and [30.417, 80]
"""
import dpg_improved as dp
import dpg_close as dc
import improved as im
from improved import mp, LN

EXT = {"EXT": dc.margin_ext}
pieces = [("26.246", "31", "0.01"), ("31", "40", "0.1"), ("40", "80", "1")]
for a, b, s in pieces:
    print(f"with Prop 5.1  [{a}, {b}] step {s}: uncertified {dp.scan(mp.mpf(a), mp.mpf(b), mp.mpf(s), EXT)}", flush=True)
for a, b, s in [("26.246", "30.001", "0.001"), ("30.417", "31", "0.001"), ("31", "40", "0.1"), ("40", "80", "1")]:
    print(f"published only [{a}, {b}] step {s}: uncertified {dp.scan(mp.mpf(a), mp.mpf(b), mp.mpf(s))}", flush=True)
print("published only, just outside the windows (expected: uncertified):",
      dp.good(mp.mpf("26.245"), mp.mpf("26.246")), dp.good(mp.mpf("30.002"), mp.mpf("30.003")),
      dp.good(mp.mpf("30.415"), mp.mpf("30.416")))


def block_bound(logy, logh, inp):
    """log of (sum of gaps >= h starting in [y, 2y]) / y, Proposition 4.2, with the best valid input."""
    ly2 = logy + LN(mp.mpf("2.02"))
    logT = LN(mp.mpf("18.6") * 100 / mp.mpf("0.005")) + logy + LN(ly2) + LN(LN(ly2)) - logh
    if logT < LN(mp.mpf("3e12")):
        return None                      # all zeros below height 3e12 are on the line; negligible
    if logT <= im.LOGT_MAX_B:
        zs = im.log_zero_sum(logy, logT, dp.c0, "B24")
    elif inp == "EXT":
        zs = im.log_zero_sum(logy, logT, dp.c0, "B24", sig=dc.SIG_EXT, cB=dc.C_EXT)
    elif inp == "B23":
        zs = dp.log_zero_sum_B23(logy, logT)
    else:
        zs = im.log_zero_sum(logy, logT, dp.c0, "RAM")
    cL = 1 / ((1 - mp.mpf("2.01") / 100) * (mp.mpf("0.999") - mp.mpf("0.005")) ** 2)
    return LN(cL * 4 * dp.K * (logT + 1)) + zs


print("top block is the largest (DPG, h = P/(2x^theta)):")
for t in ("26.246", "28", "30", "30.2", "30.417", "31", "35", "50", "80"):
    t = mp.mpf(t)
    th = dp.good(t, t, EXT)
    inp = sorted(dp.dpg_inputs(t, th, EXT))[0]
    logx = mp.e ** t
    logh = logx * (1 - th) - LN(2)
    top = block_bound(logx - LN(2), logh, inp)
    lo = logh + 2 * LN(logx)
    vals = [block_bound(lo + (logx - LN(2) - lo) * mp.mpf(i) / 400, logh, inp) for i in range(401)]
    worst = max(v for v in vals if v is not None)
    print(f"  t={mp.nstr(t, 6)} theta={th} input={inp}: top={mp.nstr(top, 6)} max over blocks={mp.nstr(worst, 6)}",
          worst <= top + abs(top) * mp.mpf("1e-12"))
