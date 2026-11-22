"""Section 9: the necessary condition of Theorem 1.6 and the one-zero barrier.

(1) Necessary condition: graphic => 2 S(g) <= P + K(g)^2 for every g.  Test on toy sequences.
(2) One-zero barrier: any L^2 / zero-density argument must at least beat the contribution of ONE
    hypothetical zero at the edge of the best known zero-free region, with the most optimistic
    pair factor (log T / 2pi) and target fraction 1/2.  Find the smallest t where that is possible.
"""
import mpmath as mp, random
mp.mp.dps = 30
LN = mp.log

def eg_graphic(d):
    d = sorted(d, reverse=True); n = len(d)
    if sum(d) % 2: return False
    for k in range(1, n + 1):
        if sum(d[:k]) > k * (k - 1) + sum(min(k, x) for x in d[k:]): return False
    return True

def necessary(d):
    P = sum(d)
    for g in set(d):
        K = sum(1 for x in d if x >= g); S = sum(x for x in d if x >= g)
        if 2 * S > P + K * K: return False
    return True

# (1) necessary condition never rejects a graphic sequence
random.seed(1); bad = 0
for _ in range(20000):
    n = random.randint(2, 14); d = [random.randint(1, n - 1) for _ in range(n)]
    if sum(d) % 2: d[0] += 1 if d[0] < n - 1 else -1
    if eg_graphic(d) and not necessary(d): bad += 1
print("graphic but violating necessary condition:", bad)
# toy 'prime-gap-like' sequence with K ~ 0.5 sqrt(P) gaps of size ~ 1.5 sqrt(P) covering > P/2
m = 400; small = [2] * (m * m // 2); K = m // 2; big = [3 * m] * K
d = small + big; P = sum(d)
print("toy: P=%d K=%d S=%d frac=%.3f necessary=%s" % (P, K, sum(big), sum(big) / P, necessary(d)))

# (2) one-zero barrier
def eta(logT, t):
    cl = 1 / (mp.mpf("5.558691") * logT)                       # Mossinghoff-Trudgian-Yang
    lw = LN(logT) / (mp.mpf("21.233") * logT)                  # Yang (Littlewood type)
    vk = 1 / (mp.mpf("53.989") * logT ** (mp.mpf(2) / 3) * LN(logT) ** (mp.mpf(1) / 3))
    return max(cl, lw, vk)

def one_zero_margin(t):
    logx = mp.e ** t; logT = logx / 2
    return -2 * eta(logT, t) * logx + LN(logT / (2 * mp.pi)) - LN(mp.mpf("0.5"))

for t in [3.72, 5, 8, 10, 13, 15, 18, 20, 22, 23]:
    print("t=%5.2f  log(one-zero bound / target) = %s" % (t, mp.nstr(one_zero_margin(mp.mpf(t)), 6)))
a, b = mp.mpf(5), mp.mpf(30)
while b - a > 1e-4:
    m_ = (a + b) / 2; (a, b) = (a, m_) if one_zero_margin(m_) < 0 else (m_, b)
print("one-zero barrier: L^2 method impossible for t <", mp.nstr(b, 6))
