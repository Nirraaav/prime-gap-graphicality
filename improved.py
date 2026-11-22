"""Graphicality threshold (Theorem 1.2): Proposition 4.2 (block estimate) + Lemma 6.1 (zero sums).

Everything is computed in log space with mpmath, since x = exp(exp(t)) with t around 25-35.
A "margin" is log(bound) - log(target); a negative margin means the inequality holds.

Block estimate (Proposition 4.2): for gaps of length >= h starting in [y, 2y],
  delta = h/(kappa y),  T = 18.6 kappa y log(2.02y) loglog(2.02y) / (eps h),
  sum of gaps <= 4 c_L K_Psi (log T + 1) S_y / y,  c_L = 1/((1 - 2.01/kappa)(0.999 - eps)^2),
where K_Psi = sum_j sup_{|tau-2j|<=1} min(F0, V/tau^2) comes from the trapezoid weight (Lemma 4.1).
The paper uses a = 0.2, kappa = 100, eps = 0.005.
"""
import mpmath as mp
import sys

mp.mp.dps = 40
LN = mp.log

# ---------------- inputs (all cited) ----------------
C0_PUBLISHED = mp.mpf("53.989")      # Bellotti JMAA 2024 VK constant (cited by Bellotti-Trudgian-Yang 2026)
# Lemma 2.6 (Bellotti 2024, App. A last range, exp(481958) < T <= exp(6.7e12)); constant recomputed with the
# corrections listed in bellotti_extend.py.  Below log T = 481958 we use Ramare only.
C_B, B_B, SIG_B = mp.mpf("6.5e11"), mp.mpf("1.448"), mp.mpf("0.9927")
LOGT_MIN_B = mp.mpf(481958)
# Lemma 3.2 (our large-sieve bound for prime pairs, see sieve_constant.py): for log x >= e^24 and even
# 2 <= h <= 80 log x, #{p <= x: p + h prime} <= SIEVE * Sg(h) x / log^2 x.
SIEVE = mp.mpf("8.001")
LOGT_MAX_B = mp.mpf("6.7e12")


def kernel_constant(a=mp.mpf("0.2"), L=LN(mp.mpf("2.01")), J=4000):
    """sum_j m_j for the trapezoid weight; sup over sigma' in [1,3] done on a grid + margin."""
    def F0():
        # int sup_{s in[1,3]} e^{s z} psi(z) dz
        left = mp.quad(lambda z: mp.e ** z * (z + a) / a, [-a, 0])
        mid = mp.quad(lambda z: mp.e ** (3 * z), [0, L])
        right = mp.quad(lambda z: mp.e ** (3 * z) * (L + a - z) / a, [L, L + a])
        return left + mid + right

    def TV(s):
        # f = e^{s z} psi(z); f' piecewise; TV(f') = jumps + int |f''|
        e = lambda z: mp.e ** (s * z)
        jumps = e(-a) / a + e(0) / a + e(L) / a + e(L + a) / a
        f2_left = mp.quad(lambda z: abs(e(z) * (s * s * (z + a) / a + 2 * s / a)), [-a, 0])
        f2_mid = mp.quad(lambda z: e(z) * s * s, [0, L])
        f2_right = mp.quad(lambda z: abs(e(z) * (s * s * (L + a - z) / a - 2 * s / a)), [L, L + a])
        return jumps + f2_left + f2_mid + f2_right

    f0 = F0()
    V = max(TV(mp.mpf(1) + mp.mpf(i) / 50) for i in range(101)) * mp.mpf("1.01")
    k = lambda tau: f0 if tau <= 0 else min(f0, V / tau ** 2)
    tot = k(0)  # j = 0 window [-1,1]
    for j in range(1, J + 1):
        tot += 2 * k(mp.mpf(2 * j - 1))
    tot += 2 * V / (2 * (2 * J + 1) - 2)  # tail: sum_{j>J} V/(2j-1)^2 <= V/(2(2J+1)-2)... (generous)
    return tot, f0, V


def log_eta(logT, c0):
    return -LN(c0) - mp.mpf(2) / 3 * LN(logT) - LN(LN(logT)) / 3


def lse(v):
    m = max(v)
    return m + LN(mp.fsum([mp.e ** (x - m) for x in v]))


def log_zero_sum(logy, logT, c0, dens="B24", sig=SIG_B, cB=C_B):
    """log of S_y / y^2, S_y = sum_{beta>=1/2,0<gamma<=T} y^{2beta} (Lemma 6.1).
    dens="B24": Ramare on [0.6, sig], N(s,T) <= cB T^{1.448(1-s)} on [sig, 1-eta].
    dens="RAM": Ramare on [0.6, 1-eta]."""
    eta = mp.e ** log_eta(logT, c0)
    lNT = logT + LN(logT)
    terms = [lNT - logy,                                            # y * N(T)
             LN(2 * logy * mp.mpf("0.1")) + lNT - mp.mpf("0.8") * logy]   # sigma in [1/2,0.6]
    rate = 2 * logy - mp.mpf(8) / 3 * (LN(3) + logT) - 2 * LN(logT)
    assert rate > 0
    if dens == "B24":
        u1 = 1 - sig
        # Ramare on [0.6, SIG_B]: integrand decreasing in u -> <= value at u1 / rate
        terms.append(LN(2 * logy * 965 / rate) + mp.mpf(8) / 3 * u1 * (LN(3) + logT)
                     + (3 + 2 * u1) * LN(logT) - 2 * u1 * logy)
        terms.append(LN(2 * logy * mp.mpf("51.5") / (2 * logy)) + 2 * LN(logT) - 2 * u1 * logy)
        R = 2 * logy - B_B * logT
        terms.append(LN(2 * cB * logy / R) - eta * R)
    elif dens == "RAM":
        terms.append(LN(2 * logy * 965 / rate) + 3 * LN(logT) - eta * rate)
        terms.append(LN(mp.mpf("51.5")) + 2 * LN(logT) - 2 * eta * logy)
    return lse(terms)


def block_norm(logy, logh, c0, dens, kappa, eps, Ksum, style="new"):
    """log of (sum of gaps >= h starting in [y,2y]) / y.  Returns (value, logT)."""
    if style == "new":
        ly2 = logy + LN(mp.mpf("2.02"))
        logT = LN(mp.mpf("18.6") * kappa / eps) + logy + LN(ly2) + LN(LN(ly2)) - logh
        cL = 1 / ((1 - mp.mpf("2.01") / kappa) * (mp.mpf("0.999") - eps) ** 2)
        K = cL * 2 * Ksum * (logT + 1) / logT
    else:  # AFGW constants (corrected two-sided pair count), as in thresholds.py
        logT = LN(8) + logy + 2 * LN(logy) - logh
        K = (2 / mp.mpf("0.249")) * 4 * 153 * 2 * mp.mpf("2.0767") * (logT + LN(3)) / logT
    if logT < LN(mp.mpf("3e12")):
        return -mp.mpf(10) ** 6, logT
    if dens == "B24" and logT <= LOGT_MIN_B:
        dens = "RAM"
    zs = log_zero_sum(logy, logT, c0, dens)
    if style == "new":     # sum over all zeros of y^{2 beta} <= 2 S_y + 2 y N(T)  (Step 5 of Proposition 4.2)
        zs = lse([zs, logT + LN(logT) - logy])
    return LN(K) + LN(logT) + zs, logT


def margin_graphic(t, c0, dens="B24", kappa=40, eps=mp.mpf("0.01"), Ksum=None, style="new",
                   theta_h=mp.mpf("0.02"), lam=40):
    """Theorem 3.5: enough that S(h0) <= eps0 * P, h0 = theta_h sqrt(P), eps0 = lam/(1+SIEVE lam) - theta_h^2 - 1e-6."""
    logx = mp.e ** t
    c1 = mp.mpf(lam) / (1 + SIEVE * lam)
    epsP = c1 - theta_h ** 2 - mp.mpf("1e-6")
    logh = LN(theta_h) + logx / 2
    top, logT = block_norm(logx - LN(2), logh, c0, dens, kappa, eps, Ksum, style)
    maxgap = LN(logx / LN(2) + 1) - LN(25) - 2 * LN(logx)       # (#blocks) * Dusart max gap / x
    return lse([top, maxgap]) - LN(epsP), logT


def threshold(f, lo=20, hi=40, tol=mp.mpf("1e-4"), **kw):
    a, b = mp.mpf(lo), mp.mpf(hi)
    assert f(b, **kw)[0] < 0 and f(a, **kw)[0] >= 0
    while b - a > tol:
        m = (a + b) / 2
        (a, b) = (a, m) if f(m, **kw)[0] < 0 else (m, b)
    return b


if __name__ == "__main__":
    # Lemma 4.1 constants and the graphicality thresholds of Table 2 (rows 3 and 4).
    Ksum, f0, V = kernel_constant(a=mp.mpf("0.2"))
    print(f"kernel (a=0.2): F0={mp.nstr(f0, 6)}  V<={mp.nstr(V, 6)}  K_Psi<={mp.nstr(Ksum, 6)}")
    cL = 1 / ((1 - mp.mpf("2.01") / 100) * (mp.mpf("0.999") - mp.mpf("0.005")) ** 2)
    print(f"c_L={mp.nstr(cL, 6)}  2 c_L K_Psi={mp.nstr(2 * cL * Ksum, 6)}  sieve constant c_S={mp.nstr(SIEVE, 6)}  "
          f"eps0={mp.nstr(40 / (1 + SIEVE * 40) - mp.mpf('0.0004') - mp.mpf('1e-6'), 8)}")
    for c0 in (C0_PUBLISHED, mp.mpf("51.34")):
        old = threshold(margin_graphic, c0=c0, style="old")
        new = threshold(margin_graphic, c0=c0, Ksum=Ksum, kappa=100, eps=mp.mpf("0.005"))
        print(f"c0={c0}: AFGW block constants t*={mp.nstr(old, 6)}   Proposition 4.2 t*={mp.nstr(new, 6)}")
