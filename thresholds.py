"""Shared inputs and margin functions; used for Theorem 1.5 and for the comparison with AFGW (Table 2).

Everything is done in log space with mpmath, since x = exp(exp(t)) with t around 25-35.
A "margin" is log(bound) - log(target); a negative margin means the inequality holds at that t.

Inputs (statements as in Section 2 of the paper):
  c0      Vinogradov-Korobov constant, beta < 1 - 1/(c0 (log t)^{2/3} (log log t)^{1/3}); 53.989 (Bellotti, JMAA 2024).
  Bellotti 2024 (arXiv:2405.12545, Thm 1.1): N(s,T) <= 1.62e11 T^{1.448(1-s)}, s in [0.9927, 1], 3e12 < T <= exp(6.7e12).
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

C_B24, B_B24, SIG_B24, LOGT_MAX_B24 = mp.mpf("1.62e11"), mp.mpf("1.448"), mp.mpf("0.9927"), mp.mpf("6.7e12")
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


def log_zero_sum_sq(logy, logT, c0, dens, s1=SIG_B24, u_switch=mp.mpf("0.0004")):
    """log of [S_y / y^2]."""
    eta = mp.e ** log_eta(logT, c0)
    terms = []
    terms.append(log_NT(logT) - logy)                            # beta = 1/2 part (S2 = y N(T))
    terms.append(log_NT(logT) - mp.mpf("0.8") * logy)            # sigma in [1/2, 0.6]
    if dens == "B24":
        # Ramare on [0.6, s1] (integrand max at s1 when rate>0), Bellotti 2024 on [s1, 1-eta]
        u1 = 1 - s1
        rateR = 2 * logy - mp.mpf(8) / 3 * (LN(3) + logT) - 2 * LN(logT)
        assert rateR > 0
        terms.append(LN(2 * logy * mp.mpf("0.4")) + LN(965) + mp.mpf(8) / 3 * u1 * (LN(3) + logT)
                     + (3 + 2 * u1) * LN(logT) - 2 * u1 * logy)
        terms.append(LN(2 * logy * mp.mpf("0.4") * mp.mpf("51.5")) + 2 * LN(logT) - 2 * u1 * logy)
        R = 2 * logy - B_B24 * logT
        assert R > 0
        terms.append(LN(2 * C_B24 * logy / R) - eta * R)
    elif dens == "RAM":
        # Ramare on [0.6, 1-eta]: integrand y^{-2u} (3T)^{8u/3} (logT)^{3+2u}, u=1-s in [eta,0.4]
        R = 2 * logy - mp.mpf(8) / 3 * (LN(3) + logT) - 2 * LN(logT)
        assert R > 0
        terms.append(LN(2 * logy * 965 / R) + 3 * LN(logT) - eta * R)
        terms.append(LN(2 * logy * mp.mpf("0.4") * mp.mpf("51.5")) + 2 * LN(logT) - 2 * eta * logy)
    elif dens == "PAPER":
        # the paper's (2.2) used for all s in [0.6,1] (NOT justified, reproduced for comparison)
        R = 2 * logy - mp.mpf(8) / 3 * logT - 2 * LN(logT)
        terms.append(LN(2 * logy * mp.mpf("2.375") / R) + 3 * LN(logT) - eta * R)
    elif dens == "B23":
        # Ramare on [0.6, 1-u_switch], Bellotti 2023 (near unity) on [1-u_switch, 1-eta]
        us = u_switch
        rateR = 2 * logy - mp.mpf(8) / 3 * (LN(3) + logT) - 2 * LN(logT)
        assert rateR > 0
        terms.append(LN(2 * logy * mp.mpf("0.4")) + LN(965) + mp.mpf(8) / 3 * us * (LN(3) + logT)
                     + (3 + 2 * us) * LN(logT) - 2 * us * logy)
        terms.append(LN(2 * logy * mp.mpf("0.4") * mp.mpf("51.5")) + 2 * LN(logT) - 2 * us * logy)
        # integrand exp(f(u)), f(u) = -2u logy + B u^{3/2} logT on u in [eta, us];
        # f'(u) <= -r with r = 2 logy - 1.5 B us^{1/2} logT > 0, so int <= exp(f(eta))/r
        f = lambda u: -2 * u * logy + B_B23 * u ** mp.mpf(1.5) * logT
        r = 2 * logy - mp.mpf(1.5) * B_B23 * mp.sqrt(us) * logT
        assert r > 0
        terms.append(LN(2 * logy * C_B23 / r) + D_B23 * LN(logT) + f(eta))
    else:
        raise ValueError(dens)
    return logsumexp(terms)


# pair-sum constant: sum_{rho'} 1/(1+|gamma-gamma'|^2) <= 2*2.0767*log(3T+1) (Lemma 2.1 used on
# both sides of gamma; the paper wrote 3 log gamma, which counts only one side).
def log_pair(logT, corrected=True):
    if corrected:
        return LN(2 * mp.mpf("2.0767") * (logT + LN(3) + mp.mpf("1e-9")))
    return LN(3 * logT)


def log_large_gap_sum_over_x(logx, logh, c0, dens, corrected=True):
    """log of an upper bound for (sum of gaps >= h with p <= x) / x, using the Heath-Brown
    argument (Lemma 3.1) on dyadic blocks [y,2y], y = x/2^{j+1}.  The block bound
    Q(y)/y is increasing in y (checked numerically in check_monotone); the sum over blocks
    is <= x * Q(x/2)/(x/2), plus the max-gap terms (<= log2(x) * x/(25 log^2 x), Dusart)."""
    logy = logx - LN(2)
    logT = LN(8) + logy + 2 * LN(logy) - logh          # T = 2 delta^{-1} log^2 y, delta = h/(4y)
    pairc = mp.e ** log_pair(logT, corrected) / logT if corrected else 3
    # I <= 4*153*pairc*logT * delta^2 y S_y ; gaps <= (2/0.249) I/(delta^2 y^2)
    K = (2 / mp.mpf("0.249")) * 4 * 153 * pairc
    main = LN(K) + LN(logT) + log_zero_sum_sq(logy, logT, c0, dens)
    maxgap = LN(logx / LN(2) + 1) - LN(25) - 2 * LN(logx)
    return logsumexp([main, maxgap]), logT


# ---------------------------------------------------------------------------
# Criteria
# ---------------------------------------------------------------------------

def margin_graphic_new(t, c0, dens, theta_h=mp.mpf("0.02"), lam=40):
    """New reduction (Erdos-Gallai cases A/B + sieve + Dusart): suffices S(h0) <= eps*P,
    h0 = theta_h*sqrt(P), eps = c1 - theta_h^2 - o(1), c1 = lam/(1+8 lam) (Siebert's sieve)."""
    logx = mp.e ** t
    c1 = mp.mpf(lam) / (1 + 8 * lam)
    eps = c1 - theta_h ** 2 - mp.mpf("1e-6")
    logh = LN(theta_h) + logx / 2
    lg, logT = log_large_gap_sum_over_x(logx, logh, c0, dens)
    return lg - LN(eps), logT


def margin_graphic_paperstyle(t, c0, dens, corrected=True):
    """The paper's reduction: S_{N'}(p_n) <= (n-1)/3 with N' = sqrt(p)/(3 log p), via our
    block bound (alpha -> 1/4 limit) -- used to quantify each improvement separately."""
    logx = mp.e ** t
    logh = logx / 2 - LN(3) - LN(logx)
    lg, logT = log_large_gap_sum_over_x(logx, logh, c0, dens, corrected)
    # target (n-1)/(3x) >= 1/(3 log x)  (pi(x) >= x/log x)
    return lg - (-LN(3) - LN(logx)), logT


def margin_paper_original(t, c0=mp.mpf("51.34"), c1=mp.mpf("2.375"), alpha=mp.mpf("0.249"), A=mp.mpf(8) / 3):
    """Exactly the paper's final inequality (Section 4.2): log(3*8.7C)+5t < (8/3)alpha*eta*log p."""
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


def margin_dpg_S(t, theta, c0, dens):
    """DPG via Lemma 2.8 with delta = floor(P/(2 x^theta)): need S(delta) <= P/4 (given d <= x^theta)."""
    logx = mp.e ** t
    logh = logx * (1 - theta) - LN(2) - mp.mpf("1e-9")
    lg, logT = log_large_gap_sum_over_x(logx, logh, c0, dens)
    return lg - LN(mp.mpf("0.2499")), logT


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


def grid_ok(fun, t0, t1, n=400, **kw):
    """Check fun(t) < 0 on a fine grid of [t0, t1]; returns (ok, worst margin)."""
    worst = None
    for i in range(n + 1):
        t = t0 + (t1 - t0) * mp.mpf(i) / n
        m = fun(t, **kw)[0]
        worst = m if worst is None else max(worst, m)
    return worst < 0, worst


if __name__ == "__main__":
    # Table 2, rows 1-2: AFGW's final inequality as printed, and with a uniform density constant.
    for c0 in (mp.mpf("53.989"), mp.mpf("51.34")):
        t = threshold(margin_paper_original, c0=c0)
        r = threshold(margin_paper_original, c0=c0, c1=mp.mpf(965))
        print(f"c0={c0}: AFGW final inequality as printed: t* = {mp.nstr(t, 6)};  with Ramare's 965: t* = {mp.nstr(r, 6)}")
    assert abs(threshold(margin_paper_original) - mp.mpf("30.313")) < 0.01
