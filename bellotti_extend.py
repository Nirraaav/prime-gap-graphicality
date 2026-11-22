"""Proposition 5.1: Bellotti's log-free zero-density constant for exp(6.7e12) < T <= exp(2e13).

Bellotti, "An explicit log-free zero density estimate for the Riemann zeta-function" (arXiv:2405.12545),
ships a program that computes the constants of her Theorem 1.1:
  https://github.com/ChiaraBellotti/Log-free-zero-density-estimates-for-zeta  (file constants.py)
That repository has no license, so this script does not redistribute it. It downloads constants.py at a
pinned commit, checks its SHA-256, and runs it with new parameters:
  * t (= T0), tupp (= T1), sigma (= alpha0), u, x, v, w;
  * fix 1: fracd2 / fracd2w call contrd2(u,v,w,x) with permuted arguments -> contrd2(u,x,v,w);
  * fix 2: cond1 uses exp(1/(6 log t)); (3.10) of the paper has e^{1/(6 gamma)} with gamma > 0.447 log T;
  * fix 3: cond1 uses (u+w)^4; log^4(VW) = (v+w)^4 log^4 T.
  * fix 4: Psisquare (d3) divides by log t once too often (Lemma 3.2 gives d3 = 4(1 + 1.0061(u L + ...))/((v-u)^2 L));
  * fix 5: constinfrontint (b1) misses a factor sqrt(2) on the r^2 terms of Corollary 2.3 with r = sqrt(2)(1 - alpha);
  * fix 6: Section 4.1 of the paper assumes that, for each picked zero, at most one other picked zero lies at distance
    < delta and at most one at each distance in [n delta, (n+1) delta).  With one picked zero per box of height delta,
    there can be up to 3 and up to 4 (two on each side).  We multiply c2 by 3 and c3, c4 by 4.
Each fix makes the constant larger or leaves it unchanged.

The RH-verification step of her proof (zeros with gamma <= 0.447 log T) is replaced, for sigma >= 0.994
and log T <= 2e13, by the classical zero-free region 1/(5.558691 log gamma); the last line printed checks it.
"""
import hashlib
import io
import os
import re
import contextlib
import urllib.request
import warnings
import mpmath as mpm

warnings.filterwarnings("ignore", category=SyntaxWarning)

COMMIT = "a43c9e5bbfa7aecc666879483d61673dd0ed88b9"
URL = f"https://raw.githubusercontent.com/ChiaraBellotti/Log-free-zero-density-estimates-for-zeta/{COMMIT}/constants.py"
SHA256 = "afd670bf9de153add60ba6c46f3a0f7f284305c461518911f118e1ac9771f104"
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".cache", "bellotti_constants.py")


def bellotti_source():
    if not os.path.exists(CACHE):
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        data = urllib.request.urlopen(URL, timeout=60).read()
        with open(CACHE, "wb") as f:
            f.write(data)
    data = open(CACHE, "rb").read()
    if hashlib.sha256(data).hexdigest() != SHA256:
        raise RuntimeError(f"{CACHE} does not match the pinned SHA-256; delete it and retry")
    return data.decode()


SRC = bellotti_source()


def constants(t_expr, tupp_expr, sigma="0.985", u="0.3520926", x="0.7238738", v="0.3616651", w="0.0000662",
              fix=True):
    s = SRC
    s = re.sub(r"^t = .*$", f"t = {t_expr}", s, flags=re.M)
    s = re.sub(r"^tupp= .*$", f"tupp= {tupp_expr}", s, flags=re.M)
    for name, val in (("u", u), ("x", x), ("v", v), ("w", w), ("sigma", sigma)):
        s = re.sub(rf"^{name}=.*$", f"{name}={val}", s, flags=re.M)
    if fix:
        s = s.replace("contrd2(u,v,w,x)", "contrd2(u,x,v,w)").replace("contrd2w(u,v,w,x)", "contrd2w(u,x,v,w)")
        old = "pow(u+w,4)*pow(2*pi,0.5)*exp(fdiv(1,6*log(t)))"
        assert s.count(old) == 1
        s = s.replace(old, "pow(v+w,4)*pow(2*pi,0.5)*exp(fdiv(1,6*0.447*log(t)))")
    ns = {}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(s, "bellotti_constants", "exec"), ns)
    U, X, V, W = (mpm.mpf(z) for z in (u, x, v, w))
    if fix:
        t, r, sq2 = ns["t"], ns["r"], mpm.sqrt(2)
        d3 = ns["Psisquare"]
        ns["Psisquare"] = lambda *a: d3(*a) * mpm.log(t) ** 2                     # fix 4
        lt = mpm.log(t + r)
        ns["CN1"] = ((sq2 / 2 + 32 * sq2 * r ** 2 / (1 + 2 * r) ** 2) * lt
                     - 16 * sq2 * r ** 2 / (1 + 2 * r) ** 2) / mpm.log(t)               # fix 5
        for name, m in (("constc2", 3), ("constc3", 4), ("constc4", 4)):           # fix 6
            f = ns[name]
            ns[name] = (lambda f, m: (lambda *a: m * f(*a)))(f, m)
    out = {k: ns[k] for k in ("CN1", "CN2")}
    for k in ("FinalC", "cond1", "condition2", "constc2", "constc3", "constc4", "Boundfirstfactor"):
        out[k] = ns[k](U, X, V, W)
    out["B"] = 2 * X
    return out


if __name__ == "__main__":
    mpm.mp.dps = 30
    rows = (("exp(481958)", "exp(6.7*power(10,12))", "Lemma 2.6 range"),
            ("exp(481958)", "exp(2*power(10,13))", "Proposition 5.1 range"))
    for T0, T1, label in rows:
        r0 = constants(T0, T1, fix=False)
        r1 = constants(T0, T1)
        print(f"{label}: T0={T0}, T1={T1}:  as shipped C = {mpm.nstr(r0['FinalC'], 6)};  "
              f"with fixes 1-6 C = {mpm.nstr(r1['FinalC'], 6)}  (B = {mpm.nstr(r1['B'], 6)})")
    # consistency: C decreases as T0 grows (constants are evaluated at T0)
    for L0 in ("6.7*power(10,12)", "power(10,13)", "2*power(10,13)"):
        print(f"  T0=exp({L0}), T1=exp(2e13): C = {mpm.nstr(constants(f'exp({L0})', 'exp(2*power(10,13))')['FinalC'], 6)}")
    b_lo = constants("exp(481958)", "exp(2*power(10,13))")["CN1"]
    b_hi = constants("exp(2*power(10,13))", "exp(2*power(10,13))")["CN1"]
    print("b1 at T0 and at T1:", mpm.nstr(b_lo, 10), mpm.nstr(b_hi, 10), " relative growth", mpm.nstr(b_hi / b_lo - 1, 3))
    print("classical ZFR at 0.447 log T1: beta <", 1 - 1 / (5.558691 * mpm.log(0.447 * 2e13)), "(need < 0.994)")
