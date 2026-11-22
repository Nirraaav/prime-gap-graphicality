# prime-gap-graphicality

Computations for the paper *Improved explicit bounds for the graphicality of the prime gap sequence*
(N. Bhattad, 2026).

Let PD_n be the sequence of the first n prime gaps, starting from p_0 = 1. Without assuming the Riemann
Hypothesis, the paper proves:

| Result | Threshold (t = log log p_n) | Script |
|---|---|---|
| Theorem 1.2: PD_n is graphic | t ≥ 25.2502 | `certify_improved.py` |
| Theorem 1.3: every realization of PD_n has (p_{n+1} − p_n)/2 independent edges | t ≥ 26.246 | `certify_dpg.py` |
| Theorem 1.5: a prime in (x, x + x^{2/3}] | log log x ≥ 24.83 | `short_interval.py` |
| Proposition 5.1: Bellotti's log-free density constant for exp(6.7·10¹²) < T ≤ exp(2·10¹³), σ ≥ 0.994 | C = 1.67·10¹¹ | `bellotti_extend.py` |
| Section 9: necessary condition and one-zero barrier | t < 20.20 | `barrier.py` |

## Setup

```
python3 -m pip install -r requirements.txt
```

Tested with Python 3.12, mpmath 1.3.0, numpy 2.2.6, scipy 1.14.1 and sympy 1.13.1. All arithmetic uses
mpmath at 40 significant digits, in log space, because x = exp(exp(t)).

## Files

- `improved.py`: Lemma 4.1 (kernel constants F0, V, K_Ψ), Proposition 4.2 (block estimate) and
  Lemma 6.1 (zero sums). Running it prints the kernel constants and the thresholds in rows 3–4 of Table 2.
- `certify_improved.py`: Theorem 1.2. It prints the thresholds of the two branches (Lemma 2.6 up to
  t = 30.2262, Ramaré from t = 29.8168), the margin Φ at the left end of each branch, a grid check that Φ is
  decreasing, and a check that the top dyadic block gives the largest bound.
- `dpg_improved.py`: the DPG margin and the grid scan, with published inputs only. Usage:
  `python3 dpg_improved.py c0 t0 t1 step`.
- `dpg_close.py`: adds Proposition 5.1 to the DPG scan for 30.0 < t < 30.43.
- `certify_dpg.py`: Theorem 1.3 on [26.246, 80], with and without Proposition 5.1.
- `dpg_hole_partial.py`: the partial DPG statement for (30.001, 30.417) when the next gap is at most
  exp(6.6999·10¹²).
- `short_interval.py`: Theorem 1.5.
- `thresholds.py`: shared margin functions, the short-interval margin, and the final inequality of
  Aggarwal–Frot–Gou–Wang (rows 1–2 of Table 2).
- `bellotti_extend.py`: Proposition 5.1.
- `barrier.py`: Theorem 1.7 tests on random sequences, and the one-zero barrier.

## Bellotti's program

`bellotti_extend.py` needs `constants.py` from
[ChiaraBellotti/Log-free-zero-density-estimates-for-zeta](https://github.com/ChiaraBellotti/Log-free-zero-density-estimates-for-zeta).
That repository has no license, so it is not copied here. The script downloads the file at commit
`a43c9e5bbfa7aecc666879483d61673dd0ed88b9` into `.cache/`, checks its SHA-256, and runs it with new
parameters and three small corrections, which are listed in the script's docstring and in Section 5 of the
paper. The first run needs network access.

## How the checks work

Each final inequality is a margin, log(bound) − log(target), and the inequality holds where the margin is
negative. For a fixed density input and a fixed θ, the margins decrease in t (Lemma 7.1 of the paper). So a
check at the left end of each range proves it on the whole range. For Theorem 1.3, each grid step [t, t'] is
certified by one (θ, input) pair that succeeds at both t and t'.

## Limitations

The scripts use floating-point arithmetic at 40 digits, not interval arithmetic. Proposition 5.1 extends a
published theorem by rerunning its author's program on a new range; it is not itself published.

## License

MIT. See `LICENSE`.
