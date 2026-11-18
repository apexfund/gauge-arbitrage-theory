# Simulation Phase Playbook: What Does V&F's Detector Respond To?

Nov 17, 2026 · @V

## Purpose

This phase finds out which kinds of market behavior make V&F's arbitrage detector read positive, negative, or zero. We run the detector on fake markets where we chose the cause, then compare its reading with a prediction and with a standard tool.

Real data can't tell us why the detector fires, because the cause is unknown. A fake market has one known behavior, so any reading can be traced back to it.

We simulate three relationships between two assets:

- **Cointegrated pair:** the gap between two prices gets pulled back toward a normal level.
- **Lead-lag:** one asset reacts to news immediately and the other reacts one step later.
- **Autocorrelation:** the spread between two assets tends to continue (momentum) or reverse (mean reversion) from one step to the next.

For each relationship, the phase delivers:

- the detector's average reading across a grid of settings, with its spread over 1,000 simulated paths
- a predicted reading computed from the true market, which shows what the detector is aiming at
- the portfolio the detector actually chose (the "calm portfolio"), which explains most surprising readings

Out of scope: real data, proofs, drift, bid-ask bounce, the 21-asset replication of V&F's Figure 2, and the power test. What to do with the results is a later section, written once results exist.

## Key ideas in plain words

The detector reads positive when a certain portfolio's returns tend to continue, and negative when they tend to reverse. Everything in this doc is about which markets produce which sign.

- **Calm portfolio:** the zero-cost, unit-length long/short mix of the two assets and cash that has the lowest variance. Cash counts as an asset that earns nothing.
- **Detector:** every few steps it re-picks the calm portfolio from the last L steps of returns. It then multiplies the portfolio's return in one step by its return in the next step and averages all those products. This is our reading of V&F's algorithm; check it against their Sections 4.1 and 4.2.
- **Same-step covariance (Γ₀):** how the two assets move together within one step.
- **Lag-1 covariance (Γ₁):** how asset i's move today relates to asset j's move next step. Rows are today, columns are next step.
- **True calm portfolio (J\*):** the calm portfolio computed from the true Γ₀ instead of from noisy data.
- **Predicted reading (D\_pred):** what the detector would average if it always used J\*.
- **Normalized reading:** D divided by the calm portfolio's variance. It is the portfolio's lag-1 autocorrelation, so it sits between −1 and 1 and lets different markets share one axis.

Optional math, for those who want it (returns have zero mean in all three relationships):

$`
\Gamma_0 = \mathrm{Cov}(r_t, r_t), \qquad \Gamma_1 = E[\, r_t \, r_{t+1}^{\top} ], \qquad D_{\mathrm{pred}} = J^{*\top} \Gamma_1 J^{*}, \qquad \rho_{\mathrm{pred}} = \frac{J^{*\top} \Gamma_1 J^{*}}{J^{*\top} \Gamma_0 J^{*}}
`$

The rest of the doc compares the detector's actual average, called D̂, against D\_pred.

## Before you start

You need a working detector function, Python with numpy and pandas, and the conventions below agreed by everyone running simulations.

**Conventions used throughout**

- Two risky assets plus cash. Cash earns zero return.
- Inputs to the detector are per-step returns (changes in log price), with zero drift. Drift is a separate test.
- The common shock has variance 1, so all other variances are relative to it.
- Every simulation is a pair of functions: `simulate` makes the returns, and `truth` returns the exact Γ₀ and Γ₁ from the settings.

**Starting defaults.** These are our choices, not V&F's. Check them against their Sections 4.1 and 4.2 before the full runs.

| Setting | Default | Meaning |
| --- | --- | --- |
| T | 5,000 | steps scored per path |
| L | 250 | window used to estimate the calm portfolio |
| Refresh | every 5 steps | how often the calm portfolio is re-picked |
| Calm directions | 1 | number of portfolios traded |
| Cash | included | cash can be part of the mix |
| N | 1,000 | paths per setting |
| Warm-up | L steps | extra steps simulated and discarded |

**Checklist before Part 1**

- [ ] The reconstructed detector runs on a (steps × 2) array of returns and returns per-step products and the weights it chose
- [ ] Everyone agrees whether the detector uses returns of prices or of log prices, and that is written in the code
- [ ] A `results/` folder exists, with one subfolder per relationship
- [ ] A fixed list of random seeds is agreed, so results can be reproduced
- [ ] Someone other than the author has read the detector against V&F Section 4.1

## Part 1: Shared infrastructure

Build four small pieces once, and all three relationships reuse them. The code below is a skeleton: the detector wrapper is the only part that depends on your own code.

**Piece 1: Market generator plus truth.** Each relationship gets a file with two functions. `simulate` makes the returns, and `truth` returns the exact Γ₀ and Γ₁ from the same settings. Check it once by simulating one path of 400,000 steps and confirming the sample covariances match `truth` to within about 0.01. This tests the formulas and the generator together.

```python
# market_<name>.py, one per relationship
import numpy as np

def simulate(p: dict, T: int, seed: int) -> np.ndarray:
    """Return a (T, 2) array of per-step returns for the two assets."""
    ...

def truth(p: dict):
    """Return (G0, G1): 2x2 same-step and lag-1 covariance matrices."""
    ...
```

**Piece 2: Prediction calculator.** From Γ₀ and Γ₁ it finds the true calm portfolio J\* and the predicted reading. Zero-cost portfolios form a flat 2-dimensional space, so the calm one is the smallest eigenvector of a 2×2 matrix. No search is needed.

```python
# predict.py
import numpy as np

# Orthonormal basis of zero-cost portfolios over (asset 1, asset 2, cash)
E = np.column_stack([np.array([1, -1, 0]) / np.sqrt(2),
                     np.array([1,  1, -2]) / np.sqrt(6)])

def predict(G0, G1):
    G = np.zeros((3, 3)); G[:2, :2] = G0          # cash has zero variance
    vals, vecs = np.linalg.eigh(E.T @ G @ E)      # eigenvalues ascending
    J = E @ vecs[:, 0]                            # lowest-variance unit mix
    var = vals[0]
    D = J[:2] @ G1 @ J[:2]                        # cash has zero return
    return {"J": J, "var": var, "D_pred": D, "rho_pred": D / var}
```

**Piece 3: Detector wrapper.** It wraps your detector so every run returns the same fields, whatever your internals look like. It must also pass a look-ahead test: the weights chosen at time t may use only returns before t.

```python
# detector_wrapper.py
def run_detector(R, L=250, refresh=5, n_dirs=1, use_cash=True) -> dict:
    """Wrap YOUR detector. Must return:
       products : per-step products of consecutive calm-portfolio returns
       port_ret : the calm portfolio's return at each step
       weights  : array (n_refresh, 3), the weights picked at each refresh
    """
    ...

def test_no_lookahead(R, t, **kw):
    """Replace every return from step t on with noise.
    Weights picked at refreshes up to t must not change."""
    rng = np.random.default_rng(0)
    R2 = R.copy(); R2[t:] = rng.standard_normal(R2[t:].shape)
    a, b = run_detector(R, **kw), run_detector(R2, **kw)
    k = ...   # number of refreshes completed at or before step t
    assert np.allclose(a["weights"][:k], b["weights"][:k])
```

**Piece 4: Runner.** It loops over paths for one setting and saves one row per path. See Part 6 for what each row holds.

## Part 2: Relationship A, cointegrated pair

When the gap between two prices is pulled back toward a normal level, we expect the detector to read negative. Its size should be −½(1−φ)²γ, where φ is the share of the gap kept each step and γ is the gap's long-run variance.

**What the market is.** Both prices share one random walk X. Asset 1 sits half the gap above it and asset 2 half the gap below. The gap z = p₁ − p₂ drifts away and is pulled back at a speed set by φ. We use this symmetric version so the two assets have equal variance, which makes the pure spread exactly the calm portfolio whenever the gap is small.

```python
# market_coint.py
def simulate(p, T, seed):
    rng = np.random.default_rng(seed)
    phi = p["phi"]
    eps = rng.normal(0, 1, T + 1)                      # shared shock, variance 1
    sd_eta = np.sqrt(p["gamma"] * (1 - phi**2)) if phi < 1 else np.sqrt(p["sig_eta2"])
    eta = rng.normal(0, sd_eta, T + 1)
    z = np.zeros(T + 1)
    z[0] = rng.normal(0, np.sqrt(p["gamma"])) if phi < 1 else 0.0
    for t in range(1, T + 1):
        z[t] = phi * z[t - 1] + eta[t]                 # the gap
    dz = np.diff(z)
    return np.column_stack([eps[1:] + dz / 2, eps[1:] - dz / 2])

def truth(p):
    phi, g = p["phi"], p.get("gamma", 0)
    if phi < 1: k, c = 2 * (1 - phi) * g, (1 - phi) ** 2 * g
    else:       k, c = p["sig_eta2"], 0.0              # random-walk gap: no reversion
    G0 = np.array([[1 + k / 4, 1 - k / 4], [1 - k / 4, 1 + k / 4]])
    G1 = (-c / 4) * np.array([[1, -1], [-1, 1]])
    return G0, G1
```

**The prediction.** When the spread is the calm portfolio, the reading and its normalized form are:

$`
D_{\mathrm{pred}} = -\tfrac{1}{2}(1-\varphi)^2 \gamma, \qquad \rho_{\mathrm{pred}} = -\frac{1-\varphi}{2}, \qquad \text{valid while } (1-\varphi)\gamma < \tfrac{2}{3}
`$

Above that limit the calm portfolio switches to both assets against cash, and the reading goes to about zero.

**Settings to run**

| Sweep | Fixed | Varied |
| --- | --- | --- |
| Reversion speed | γ = 0.25 | φ = 0, 0.5, 0.8, 0.9, 0.95, 0.99 |
| Gap size | φ = 0.9 | γ = 0.05, 0.25, 1, 5, 20 |
| Null (no reversion) | σ\_η² = 0.05 | φ = 1 |

Expected values for the reversion-speed sweep (γ = 0.25):

| φ | D\_pred | Normalized |
| --- | --- | --- |
| 0 | −0.125 | −0.50 |
| 0.5 | −0.0313 | −0.25 |
| 0.8 | −0.005 | −0.10 |
| 0.9 | −0.00125 | −0.05 |
| 0.95 | −0.00031 | −0.025 |
| 0.99 | −0.000013 | −0.005 |

In the gap-size sweep, γ = 5 still follows the formula (D\_pred = −0.005·γ), and the formula stops applying between γ = 5 and γ = 20, at γ = 6.67. The φ = 0.99 row is the hardest to see: expect it to sit only about ten standard errors from zero. The null row should read zero.
