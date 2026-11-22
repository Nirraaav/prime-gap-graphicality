"""Certification of Theorem 1.2 (graphicality for log log p_n >= 25.2502).

Usage:  python3 certify_improved.py [c0]        (default c0 = 53.989)

1. Thresholds: Lemma 2.6 branch (valid while log T <= 6.7e12) and Ramare branch (no upper limit).
2. Monotonicity: Phi(t) is nonincreasing on a grid of each branch (Lemma 7.1 proves this; this is a check).
   With monotonicity, Phi < 0 at the left end of a branch gives Phi < 0 on the whole branch.
3. Top block: for several t, the block bound B(y)/y is largest at the top block y = x/2, scanning
   y from h0 log^2 x (below Dusart's cutoff ~ 3 h0 log^2 x) up to x/2.
"""
import sys
from improved import *

c0 = mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else C0_PUBLISHED
K = kernel_constant(a=mp.mpf("0.2"))[0]
kw = dict(c0=c0, Ksum=K, kappa=100, eps=mp.mpf("0.005"))
print("target eps0 = 40/321 - 0.0004 - 1e-6 =", mp.nstr(mp.mpf(40) / 321 - mp.mpf("0.0004") - mp.mpf("1e-6"), 10))

tB = threshold(margin_graphic, dens="B24", **kw)
tR = threshold(margin_graphic, dens="RAM", **kw)
a, b = mp.mpf(25), mp.mpf(35)
while b - a > mp.mpf("1e-4"):
    m = (a + b) / 2
    (a, b) = (m, b) if margin_graphic(m, dens="B24", **kw)[1] <= LOGT_MAX_B else (a, m)
tcap = a
print(f"c0={c0}: Lemma 2.6 branch t*={mp.nstr(tB, 7)} (valid to t={mp.nstr(tcap, 7)});  Ramare branch t*={mp.nstr(tR, 7)}")
print("  Phi(25.2502) =", mp.nstr(margin_graphic(mp.mpf("25.2502"), dens="B24", **kw)[0], 6))
print("  Phi(29.8168) [Ramare] =", mp.nstr(margin_graphic(mp.mpf("29.8168"), dens="RAM", **kw)[0], 6))


def mono(t0, t1, n, dens):
    vals = [margin_graphic(t0 + (t1 - t0) * mp.mpf(i) / n, dens=dens, **kw)[0] for i in range(n + 1)]
    return max(vals[i + 1] - vals[i] for i in range(n)), vals[0], max(vals)


for (t0, t1, n, d) in ((mp.mpf("25.2502"), tcap, 800, "B24"), (mp.mpf("29.8168"), mp.mpf(80), 2000, "RAM")):
    inc, first, worst = mono(t0, t1, n, d)
    print(f"  {d} grid [{mp.nstr(t0, 6)}, {mp.nstr(t1, 6)}]: Phi(left)={mp.nstr(first, 5)}  worst={mp.nstr(worst, 5)}  "
          f"largest step change={mp.nstr(inc, 5)} (<= 0 means nonincreasing)")
for t in (100, 150, 200, 300):
    print(f"  Phi({t}) [Ramare] =", mp.nstr(margin_graphic(mp.mpf(t), dens="RAM", **kw)[0], 6))


def topmax(t, dens, n=400):
    logx = mp.e ** t
    logh = LN(mp.mpf("0.02")) + logx / 2
    top = block_norm(logx - LN(2), logh, c0, dens, 100, mp.mpf("0.005"), K)[0]
    lo = logh + 2 * LN(logx)
    worst = None
    for i in range(n + 1):
        ly = lo + (logx - LN(2) - lo) * mp.mpf(i) / n
        try:
            v = block_norm(ly, logh, c0, dens, 100, mp.mpf("0.005"), K)[0]
        except (AssertionError, ZeroDivisionError, ValueError):
            continue
        worst = v if worst is None else max(worst, v)
    return worst <= top * (1 - mp.mpf("1e-12")) if top < 0 else worst <= top


print("  top block is the largest:", [topmax(t, "B24") for t in (tB, (tB + tcap) / 2, tcap)],
      [topmax(t, "RAM") for t in (tR, mp.mpf(40), mp.mpf(80))])
