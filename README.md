# prime-gap-graphicality

Computations for the paper *Improved explicit bounds for the graphicality of the prime gap sequence*
(N. Bhattad, 2026).

Let PD_n be the sequence of the first n prime gaps, starting from p_0 = 1. Without assuming the Riemann
Hypothesis, the paper proves:

| Result | Threshold | Script |
|---|---|---|
| Theorem 1.2: PD_n is graphic | log log p_n ≥ 25.2894 | `certify_improved.py` |
| Theorem 1.3: every realization of PD_n has (p_{n+1} − p_n)/2 independent edges | log log p_n ≥ 26.285 | `find_dpg_thresholds.py`, `certify_dpg.py` |
| Theorem 1.5: a prime in (x, x + x^{2/3}] | log log x ≥ 24.978 | `short_interval.py` |
| Lemma 3.2: #{p ≤ x : p + h prime} ≤ c_S 𝔖(h) x / log² x for log x ≥ e^24 | c_S = 8.001 | `sieve_constant.py` |
| Lemma 2.6: Bellotti's log-free density constant, corrected | C = 6.5·10¹¹ | `bellotti_extend.py` |
| Proposition 5.1: the same estimate for exp(6.7·10¹²) < T ≤ exp(2·10¹³), σ ≥ 0.994 | C = 6.7·10¹¹ | `bellotti_extend.py` |
| Section 9: necessary condition and one-zero obstruction | t < 20.20 | `barrier.py` |

Theorem 1.4 is proved by hand from the output of `certify_improved.py`.

## Setup

```
python3 -m pip install -r requirements.txt
```

Tested with Python 3.12, mpmath 1.3.0, numpy 2.2.6, scipy 1.14.1 and sympy 1.13.1. All arithmetic uses
mpmath at 40 significant digits, in log space, because x = exp(exp(t)). `certify_dpg.py` takes a few minutes;
the other scripts take seconds to a minute.

## Files

- `improved.py`: Lemma 4.1 (kernel constants F0, V0, K_Ψ), Proposition 4.2 (block estimate) and Lemma 6.1
  (zero sums). Running it prints the kernel constants, ε0, and rows 3–4 of Table 2.
- `certify_improved.py`: Theorem 1.2. It prints the thresholds of the two branches (Lemma 2.6 up to
  t = 30.2262, Ramaré from t = 29.7988), the margin Φ at the left end of each branch, a grid check that Φ does not
  increase, and a check that the top dyadic block gives the largest bound.
- `dpg_improved.py`: the DPG margin, the short-interval check for each θ, and the step-by-step scan, with
  published inputs only. Usage: `python3 dpg_improved.py c0 t0 t1 step`.
- `dpg_close.py`: adds Proposition 5.1 to the DPG scan.
- `find_dpg_thresholds.py`: locates 26.285, 30.001 and 30.412 by bisection.
- `certify_dpg.py`: Theorem 1.3 on [26.285, 80], with and without Proposition 5.1, plus the top-block check.
- `short_interval.py`: Theorem 1.5.
- `sieve_constant.py`: the constant c_S of Lemma 3.2.
- `thresholds.py`: the short-interval margin and the final inequality of Aggarwal–Frot–Gou–Wang (rows 1–2 of
  Table 2).
- `bellotti_extend.py`: Lemma 2.6 (corrected constant) and Proposition 5.1.
- `barrier.py`: tests of Theorem 1.6 on random sequences, and the one-zero obstruction.

## Bellotti's program

`bellotti_extend.py` needs `constants.py` from
[ChiaraBellotti/Log-free-zero-density-estimates-for-zeta](https://github.com/ChiaraBellotti/Log-free-zero-density-estimates-for-zeta).
That repository has no license, so it is not copied here. The script downloads the file at commit
`a43c9e5bbfa7aecc666879483d61673dd0ed88b9` into `.cache/`, checks its SHA-256, and runs it with the corrections
listed in its docstring and in Section 5 of the paper. The first run needs network access.

## How the checks work

Each final inequality is a margin, log(bound) − log(target), and the inequality holds where the margin is
negative. For a fixed density input and a fixed θ, a negative margin stays negative as t grows (Lemma 7.1 of the
paper). So a check at the left end of each range proves it on the whole range. For Theorem 1.3, each grid step
[t, t'] is certified by one (θ, input) pair that succeeds at both t and t'.

## Limitations

The scripts use floating-point arithmetic at 40 digits, not interval arithmetic. The corrected constant in
Lemma 2.6 and Proposition 5.1 come from rerunning Bellotti's program with corrections; they are not published
elsewhere.

## License

MIT. See `LICENSE`.
