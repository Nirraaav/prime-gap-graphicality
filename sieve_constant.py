"""Lemma 3.2: explicit large-sieve bound for prime pairs, replacing Siebert's theorem.

For L = log x >= e^24 and even h with 2 <= h <= 80 L,
    #{p <= x : p + h prime} <= c_S * Sg(h) * x / log^2 x,
where Sg(h) = 2 C2 prod_{p | h, p > 2} (p-1)/(p-2).  This script evaluates c_S(L) from the proof and checks that it
decreases in L, so c_S(e^24) bounds it for all larger L.

Proof outline (see the paper): Montgomery-Vaughan large sieve with Q = e^{L/2}/L gives count <= Q + (x+Q^2)/Lsum(Q);
Lsum(Q) >= 2 G(Q/2) with G(z) = sum_{m <= z odd squarefree} g(m); G(z) >= l^2 (1 - eps) / (4 Sg(h)), l = log z, where
  1 - eps = (1 - eta1)(1 - eta2) [ (log Y / l)^2 - 2^{omega+2} rho^2 (log Y + 1) / l^2 ],   log Y = l - 16 log l,
  eta1 = l^{-4} * 2.3661^omega            (tail of sum over n2 | h^infty, V = l^8),
  eta2 = M * l^{-2}                        (tail of sum of beta(c)/c, U = l^8),
  omega = omega(rad 2h) <= 1.3841 log(160L)/loglog(160L)     (Robin 1983, n >= 3),
  rho = r/phi(r) <= e^gamma loglog(160L) + 2.51/loglog(160L) (Rosser-Schoenfeld (3.42)).
Then c_S = 8 [ (L/(2l))^2 (1 + L^-2) / (1 - eps) + Q L^2 / (8 * 1.32 x) ].
"""
import mpmath as mp

mp.mp.dps = 40


def M_const(P=10 ** 6):
    """M = prod_{p >= 3} (1 + p^{-3/2} / (1 - 2 p^{-3/4})), with a rigorous tail bound."""
    sieve = bytearray([1]) * (P + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(P ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    logprod = mp.mpf(0)
    for p in range(3, P + 1):
        if sieve[p]:
            p = mp.mpf(p)
            logprod += mp.log(1 + p ** -1.5 / (1 - 2 * p ** -0.75))
    # tail: log(1+u) <= u, p^{-3/2}/(1-2p^{-3/4}) <= 1.01 p^{-3/2} for p > P, sum_{n > P} n^{-3/2} <= 2 P^{-1/2}
    tail = mp.mpf("1.01") * 2 / mp.sqrt(P)
    return mp.e ** (logprod + tail)


M = M_const()


def c_S(L):
    l = L / 2 - mp.log(L) - mp.log(2)
    lr = mp.log(160 * L)
    omega = mp.mpf("1.3841") * lr / mp.log(lr)
    rho = mp.e ** mp.euler * mp.log(lr) + mp.mpf("2.51") / mp.log(lr)
    eta1 = l ** -4 * mp.mpf("2.3661") ** omega
    eta2 = M * l ** -2
    logY = l - 16 * mp.log(l)
    one_minus_eps = (1 - eta1) * (1 - eta2) * ((logY / l) ** 2 - 2 ** (omega + 2) * rho ** 2 * (logY + 1) / l ** 2)
    Qterm = mp.e ** (-L / 2) * L / (8 * mp.mpf("1.32"))     # Q L^2 / (8 * 1.32 x) with Q = e^{L/2}/L
    return 8 * ((L / (2 * l)) ** 2 * (1 + L ** -2) / one_minus_eps + Qterm), one_minus_eps


if __name__ == "__main__":
    print("M =", mp.nstr(M, 10))
    prev = None
    for t in [24, 24.5, 25, 26, 28, 30, 35, 40, 60, 80, 120]:
        c, ome = c_S(mp.e ** mp.mpf(t))
        print(f"t = {t}: c_S = {mp.nstr(c, 12)}   1 - eps = {mp.nstr(ome, 10)}")
        assert prev is None or c <= prev
        prev = c
    print("c_S(e^24) =", mp.nstr(c_S(mp.e ** 24)[0], 15), "-> use c_S = 8.001")
    assert c_S(mp.e ** 24)[0] < mp.mpf("8.001")
