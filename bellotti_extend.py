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
    out = {k: ns[k] for k in ("CN1", "CN2")}
    for k in ("FinalC", "cond1", "condition2", "constc2", "constc3", "constc4", "Boundfirstfactor"):
        out[k] = ns[k](U, X, V, W)
    out["B"] = 2 * X
    return out


if __name__ == "__main__":
    mpm.mp.dps = 30
    # 1. reproduce her stated last row: T in (exp(481958), exp(6.7e12)], alpha0 = 0.985 -> C = 1.62e11, B = 1.448
    r = constants("exp(481958)", "exp(6.7*power(10,12))", fix=False)
    print("reproduce (unfixed):  C =", mpm.nstr(r["FinalC"], 6), " B =", mpm.nstr(r["B"], 6))
    r = constants("exp(481958)", "exp(6.7*power(10,12))")
    print("reproduce (fixed):    C =", mpm.nstr(r["FinalC"], 6))
    # 2. new range: T0 = exp(6.7e12), T1 = exp(2e13)
    for L0 in ("6.7*power(10,12)", "power(10,13)", "1.5*power(10,13)", "2*power(10,13)"):
        r = constants(f"exp({L0})", "exp(2*power(10,13))")
        print(f"T0=exp({L0}), T1=exp(2e13): C = {mpm.nstr(r['FinalC'],6)}  c1 = {mpm.nstr(r['cond1'],8)}  "
              f"c5 = {mpm.nstr(r['condition2'],6)}  c3 = {mpm.nstr(r['constc3'],6)}")
    # 3. same T0 = exp(481958) with T1 = exp(2e13): one constant for the whole range
    r = constants("exp(481958)", "exp(2*power(10,13))")
    print("T in (exp(481958), exp(2e13)]: C =", mpm.nstr(r["FinalC"], 6))
    # 4. low-zero step: zeros with gamma <= 0.447 log T1 have beta < 1 - 1/(5.558691 log(0.447 log T1))
    L1 = 2e13
    print("classical ZFR at 0.447 log T1: beta <", 1 - 1 / (5.558691 * mpm.log(0.447 * L1)), "(need < 0.994)")
    # 5. b1 (CN1) is the one constant that grows with T; compare its value at T0 = exp(481958) and at T1
    b_lo = constants("exp(481958)", "exp(2*power(10,13))")["CN1"]
    b_hi = constants("exp(2*power(10,13))", "exp(2*power(10,13))")["CN1"]
    print("b1 at T0 and at T1:", mpm.nstr(b_lo, 12), mpm.nstr(b_hi, 12), " relative growth", mpm.nstr(b_hi / b_lo - 1, 3))
