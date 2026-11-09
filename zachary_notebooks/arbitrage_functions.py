"""Simulation and arbitrage detection functions developed in 11_02_2026.ipynb."""

import numpy as np
import pandas as pd


def brownian_increments(n_steps, n_shocks, dt, rng):
    """Generate independent Brownian increments with shape (steps, shocks)."""
    # Use positive integer step/shock counts, dt > 0, and a NumPy random generator.
    return rng.normal(0.0, np.sqrt(dt), size=(n_steps, n_shocks))


def brownian_paths(increments):
    """Convert a (steps, shocks) array into paths, including the initial zero."""
    # Use finite increments shaped (n_steps, n_shocks).
    increments = np.asarray(increments, dtype=float)
    initial = np.zeros((1, increments.shape[1]))
    return np.vstack([initial, np.cumsum(increments, axis=0)])


def simulate_gbm(initial_prices, alpha, sigma, dt, increments, assets=None):
    """Simulate constant-coefficient GBM; sigma has shape (assets, shocks)."""
    # Use finite initial_prices and alpha vectors of length N; initial prices and dt > 0.
    # sigma: (N, n_shocks); increments: (n_steps, n_shocks).
    # If supplied, assets should contain N column names.
    initial_prices = np.asarray(initial_prices, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    increments = np.asarray(increments, dtype=float)

    variance = np.sum(sigma**2, axis=1)
    log_r = (alpha - 0.5 * variance) * dt + increments @ sigma.T
    log_growth = np.vstack([np.zeros(len(initial_prices)), np.cumsum(log_r, axis=0)])
    return pd.DataFrame(
        initial_prices * np.exp(log_growth),
        index=pd.RangeIndex(len(increments) + 1, name="Step"),
        columns=assets,
    )


def calculate_returns(prices):
    """Return log and simple returns, labeled by their ending observation."""
    # Use a DataFrame of finite positive prices: time rows, asset columns.
    # Use unique asset names and unique, increasing observation labels.
    ratios = prices / prices.shift(1)
    return np.log(ratios).iloc[1:], (ratios - 1).iloc[1:]


def estimate_covariance(window):
    """Estimate log-return covariance per step, using denominator L (Eq. 93)."""
    # Use finite log returns shaped (L, N), with at least two rows.
    values = np.asarray(window, dtype=float)
    R = values - values.mean(axis=0)
    return R.T @ R / len(R)


def project_covariance(Omega):
    """Construct G = P Omega P, equivalent to (Eq. 85)."""
    # Use a finite symmetric covariance matrix shaped (N, N), with N >= 2.
    Omega = np.asarray(Omega, dtype=float)
    N = len(Omega)
    P = np.eye(N) - np.ones((N, N)) / N
    G = P @ Omega @ P
    return (G + G.T) / 2


def estimate_directions(Omega, k):
    """Select k orthonormal zero-cost directions from G (Eqs. 86-88)."""
    # Use a covariance matrix and an integer k between 1 and N - 1.
    G = project_covariance(Omega)
    N = len(G)
    _, eigenvectors = np.linalg.eigh(G)
    directions = []
    for vector in eigenvectors.T:
        candidate = vector - vector.mean()
        for direction in directions:
            candidate -= (candidate @ direction) * direction
        length = np.linalg.norm(candidate)
        if length > 1e-8:
            directions.append(candidate / length)
        if len(directions) == k:
            break

    J = np.column_stack(directions)
    # Recompute eigenvalues after removing the trivial component.
    eigenvalues = np.diag(J.T @ G @ J)
    return J, eigenvalues


def align_directions(current, previous=None):
    """Align directions with an SVD solution of Appendix A (Eqs. 120–122)."""
    # Use orthonormal columns shaped (N, k), with 1 <= k < N.
    # previous should have the same shape, or be None for the first basis.
    current = np.asarray(current, dtype=float)
    if previous is None:
        aligned = current.copy()
        # Treat near-ties consistently so roundoff cannot choose the opposite sign.
        magnitudes = np.abs(aligned)
        largest = magnitudes.max(axis=0)
        pivots = np.argmax(np.isclose(magnitudes, largest, rtol=1e-10, atol=0), axis=0)
        signs = np.where(aligned[pivots, np.arange(aligned.shape[1])] < 0, -1, 1)
        return aligned * signs

    previous = np.asarray(previous, dtype=float)
    overlap = current.T @ previous  # Eq. 120
    # Singular value decomposition: overlap = A @ diag(singular_values) @ Bt.
    A, singular_values, Bt = np.linalg.svd(overlap)
    # Keep the rotation/reflection and discard the scaling to align the bases.
    return current @ (A @ Bt)


def run_detector(prices, L, dt, k=1, rolling=True):
    """Run the section 4.1 detector (Eqs. 85-90 and 93)."""
    # Use the price DataFrame described in calculate_returns, and dt > 0.
    # L is an integer >= 2; leave at least two returns after the L-return window.
    # Use integer 1 <= k < N for either fixed or rolling estimation.
    log_r, simple_r = calculate_returns(prices)

    previous = None
    drifts, directions, eigenvalues, variances = [], [], [], []
    for i in range(L, len(log_r)):
        if rolling or i == L:
            Omega = estimate_covariance(log_r.iloc[i - L:i])
            current, selected_values = estimate_directions(Omega, k)
            current = align_directions(current, previous)
        drifts.append(current.T @ simple_r.iloc[i].to_numpy() / dt)
        directions.append(current.copy())
        eigenvalues.append(selected_values.copy())
        variances.append(np.diag(current.T @ Omega @ current))
        previous = current

    index = log_r.index[L:]
    names = [f"Direction {a + 1}" for a in range(k)]
    drift_table = pd.DataFrame(drifts, index=index, columns=names)
    detector = (drift_table.shift(1) * drift_table).sum(axis=1, min_count=k)
    results = pd.DataFrame({"Detector": detector})
    results["Running mean"] = detector.expanding().mean()
    return {
        "results": results,
        "drifts": drift_table,
        "directions": np.asarray(directions),
        "eigenvalues": pd.DataFrame(eigenvalues, index=index, columns=names),
        "direction_variances": pd.DataFrame(variances, index=index, columns=names),
    }


def change_numeraire(prices, asset):
    """Express each asset price in units of the selected asset."""
    # Use finite positive prices in a DataFrame; asset must be one of its columns.
    return prices.div(prices[asset], axis=0)


def check_numeraires(prices, L, dt, k=1, rolling=True):
    """Run detection in every pricing unit and retain the cross-unit range."""
    # Use the same prices, L, dt, k, and rolling inputs as run_detector.
    runs = {}
    for asset in prices.columns:
        relative_prices = change_numeraire(prices, asset)
        runs[asset] = run_detector(relative_prices, L, dt, k, rolling)
    signals = pd.DataFrame({asset: run["results"]["Detector"]
                            for asset, run in runs.items()})
    summary = pd.DataFrame({
        "Mean": signals.mean(axis=1),
        "Minimum": signals.min(axis=1),
        "Maximum": signals.max(axis=1),
    })
    summary["Range"] = summary["Maximum"] - summary["Minimum"]
    return {"runs": runs, "signals": signals, "summary": summary}
