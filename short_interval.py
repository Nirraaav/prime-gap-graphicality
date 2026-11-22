"""Certification of Theorem 1.5 (a prime in (x, x + x^{2/3}] for log log x >= 24.83).

The margin is nonincreasing in t for each density input (Lemma 7.1), and each input is valid on an
interval of t. So it suffices to evaluate the margin at the left end of each input's range:
  Lemma 2.6 (Bellotti 2024):   from t = 24.83 while log T <= 6.7e12  (up to t ~ 30.63)
  Proposition 5.1 (extended):  where 6.7e12 < log T <= 2e13           (t ~ 30.63 to ~ 31.7)
  Lemma 2.7 (Bellotti 2023):   from t = 30.95, all larger t
A grid check of monotonicity is printed as well.
"""
import thresholds as th
from thresholds import mp, LN

c0, theta = mp.mpf("53.989"), mp.mpf(2) / 3
EXT = dict(sig=mp.mpf("0.994"), cB=mp.mpf("1.67e11"))


def m(t, dens, **kw):
    return th.margin_short_interval(mp.mpf(t), theta=theta, c0=c0, dens=dens, **kw)


def t_of_logT(L):
    a, b = mp.mpf(25), mp.mpf(35)
    while b - a > mp.mpf("1e-9"):
        mid = (a + b) / 2
        (a, b) = (mid, b) if m(mid, "RAM")[1] <= L else (a, mid)
    return a


t_cap, t_ext = t_of_logT(th.LOGT_MAX_B24), t_of_logT(mp.mpf("2e13"))
print(f"log T = 6.7e12 at t = {mp.nstr(t_cap, 8)};  log T = 2e13 at t = {mp.nstr(t_ext, 8)}")
print("Lemma 2.6 at t=24.83:", mp.nstr(m("24.83", "B24")[0], 6), "  at t=24.82:", mp.nstr(m("24.82", "B24")[0], 6))
print("Prop 5.1 at the left end of its range:", mp.nstr(m(t_cap + mp.mpf("1e-8"), "B24", **EXT)[0], 6))
print("Lemma 2.7 at t=30.95:", mp.nstr(m("30.95", "B23")[0], 6), "  at t=30.94:", mp.nstr(m("30.94", "B23")[0], 6))


def mono(t0, t1, n, dens, **kw):
    v = [m(t0 + (t1 - t0) * mp.mpf(i) / n, dens, **kw)[0] for i in range(n + 1)]
    return max(v[i + 1] - v[i] for i in range(n))


print("largest step change (<= 0 means nonincreasing):",
      "Lemma 2.6", mp.nstr(mono(mp.mpf("24.83"), t_cap, 500, "B24"), 4),
      " Prop 5.1", mp.nstr(mono(t_cap + mp.mpf("1e-8"), t_ext, 200, "B24", **EXT), 4),
      " Lemma 2.7", mp.nstr(mono(mp.mpf("30.95"), mp.mpf("80"), 500, "B23"), 4))
