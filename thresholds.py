"""Shared inputs and margin functions; used for Theorem 1.5 and for the comparison with AFGW (Table 2).

Everything is done in log space with mpmath, since x = exp(exp(t)) with t around 25-35.
A "margin" is log(bound) - log(target); a negative margin means the inequality holds at that t.

Inputs (statements as in Section 2 of the paper):
  c0      Vinogradov-Korobov constant, beta < 1 - 1/(c0 (log t)^{2/3} (log log t)^{1/3}); 53.989 (Bellotti, JMAA 2024).
  Bellotti 2024 (arXiv:2405.12545, Thm 1.1, App. A last range): N(s,T) <= C T^{1.448(1-s)} for s in [0.985, 1],
          exp(481958) < T <= exp(6.7e12); published C = 1.62e11, we use the corrected C = 6.5e11 (bellotti_extend.py).
  Bellotti 2023 (arXiv:2311.05136, Thm 1.3): N(s,T) <= 4.45e12 T^{57.8875 (1-s)^{3/2}} (log T)^{10393/900},
          s in [0.98, 1], T >= exp(6.7e12).
  Ramare, corrected form quoted in Kadiri-Lumley-Ng (1.5): N(s,T) <= 965 (3T)^{8(1-s)/3} (log T)^{5-2s} + 51.5 (log T)^2,
          s >= 0.52, T >= 2000.
  margin_paper_original reproduces the final inequality of AFGW, Section 4.2, which uses the KLN Table 1 constant
  2.375 (row sigma_0 = 0.60) for all sigma; see Appendix A of the paper.
"""
import mpmath as mp

mp.mp.dps = 40
LN = mp.log
E = mp.e

# Lemma 2.6: Bellotti 2024, last range of App. A, constant recomputed with the corrections in bellotti_extend.py
C_B24, B_B24, SIG_B24, LOGT_MAX_B24 = mp.mpf("6.5e11"), mp.mpf("1.448"), mp.mpf("0.9927"), mp.mpf("6.7e12")
C_B23, B_B23, D_B23 = mp.mpf("4.45e12"), mp.mpf("57.8875"), mp.mpf(10393) / 900


def log_eta(logT, c0):
    """log of eta(T) = 1/(c0 (log T)^{2/3} (log log T)^{1/3})."""
    return -LN(c0) - mp.mpf(2) / 3 * LN(logT) - LN(LN(logT)) / 3


def logsumexp(vals):
    m = max(vals)
    return m + LN(mp.fsum([mp.e ** (v - m) for v in vals]))


# ---------------------------------------------------------------------------
# Zero sums.  For a block y and height T we bound
#   S_y / y^2  with  S_y = sum_{beta>=1/2, 0<gamma<=T} y^{2 beta}
#   = y N(T) + int_{1/2}^{1} 2 log y * y^{2s} N(s,T) ds            (s = sigma)
# returning log of the bound.  `dens` selects the near-1 density estimate.
# ---------------------------------------------------------------------------

def log_NT(logT):
    # N(T) <= T log T  (crude, valid for T >= 3e12)
    return logT + LN(logT)


def margin_paper_original(t, c0=mp.mpf("51.34"), c1=mp.mpf("2.375"), alpha=mp.mpf("0.249"), A=mp.mpf(8) / 3):
    """Exactly the final inequality of AFGW (Section 4.2): log(3*8.7C)+5t < (8/3)alpha*eta*log p."""
    C = 4 * 1836 * c1 / (A * alpha) * (2 / A - alpha) ** 4
    logp = mp.e ** t
    L1 = logp / 2 + 3 * t + LN(12)
    rhs = A * alpha / c0 * logp / (L1 ** (mp.mpf(2) / 3) * LN(L1) ** (mp.mpf(1) / 3))
    return LN(3 * mp.mpf("8.7") * C) + 5 * t - rhs, L1


def margin_short_interval(t, theta, c0, dens, sig=SIG_B24, cB=C_B24):
    """Prime in (x, x + x^theta]: need 2*sum_{0<gamma<T} x^{beta-1} + 0.01 + pp/h < 1."""
    logx = mp.e ** t
    logh = theta * logx
    # T = 920 (x+h) log(x+h) loglog(x+h)/h
    logT = LN(920) + (1 - theta) * logx + LN(logx) + LN(LN(logx)) + mp.mpf("1e-9")
    eta = mp.e ** log_eta(logT, c0)
    terms = [log_NT(logT) - logx / 2,                                  # zeros with beta <= 1/2
             log_NT(logT) - mp.mpf("0.4") * logx + LN(mp.mpf("0.1") * logx)]  # [1/2,0.6]
    if dens in ("B24", "B23"):
        s1 = sig if dens == "B24" else 1 - mp.mpf("0.0004")
        u1 = 1 - s1
        rate = logx - mp.mpf(8) / 3 * (LN(3) + logT) - 2 * LN(logT)
        assert rate > 0, "Ramare middle range has no saving"
        terms.append(LN(mp.mpf("0.4") * logx * 965) + mp.mpf(8) / 3 * u1 * (LN(3) + logT)
                     + (3 + 2 * u1) * LN(logT) - u1 * logx)
        terms.append(LN(mp.mpf("0.4") * logx * mp.mpf("51.5")) + 2 * LN(logT) - u1 * logx)
        if dens == "B24":
            R = logx - B_B24 * logT
            terms.append(LN(cB * logx / R) - eta * R)
        else:
            f = lambda u: -u * logx + B_B23 * u ** mp.mpf(1.5) * logT
            r = logx - mp.mpf(1.5) * B_B23 * mp.sqrt(u1) * logT
            assert r > 0
            terms.append(LN(logx * C_B23 / r) + D_B23 * LN(logT) + f(eta))
    elif dens == "RAM":
        R = logx - mp.mpf(8) / 3 * (LN(3) + logT) - 2 * LN(logT)
        terms.append(LN(logx * 965 / R) + 3 * LN(logT) - eta * R)
    tot = LN(2) + logsumexp(terms)
    pp = LN(4) + LN(logx) + (-logx / 2)          # prime powers / h, generous
    lhs = logsumexp([tot, LN(mp.mpf("0.01")), pp])
    return lhs - LN(mp.mpf("0.98")), logT


def threshold(fun, lo=15, hi=40, tol=mp.mpf("1e-4"), **kw):
    """Smallest t in [lo,hi] (to tol) with fun(t) < 0, assuming sign change once."""
    flo = fun(mp.mpf(lo), **kw)[0]
    fhi = fun(mp.mpf(hi), **kw)[0]
    if flo < 0:
        return mp.mpf(lo)
    if fhi >= 0:
        return None
    a, b = mp.mpf(lo), mp.mpf(hi)
    while b - a > tol:
        m = (a + b) / 2
        if fun(m, **kw)[0] < 0:
            b = m
        else:
            a = m
    return b


if __name__ == "__main__":
    # Table 2, rows 1-2: AFGW's final inequality as printed, and with a uniform density constant.
    for c0 in (mp.mpf("53.989"), mp.mpf("51.34")):
        t = threshold(margin_paper_original, c0=c0)
        r = threshold(margin_paper_original, c0=c0, c1=mp.mpf(965))
        print(f"c0={c0}: AFGW final inequality as printed: t* = {mp.nstr(t, 6)};  with Ramare's 965: t* = {mp.nstr(r, 6)}")
    assert abs(threshold(margin_paper_original) - mp.mpf("30.313")) < 0.01
