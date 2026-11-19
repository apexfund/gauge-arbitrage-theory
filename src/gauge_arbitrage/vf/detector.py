import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.linalg import null_space

#Log returns (since the detector run_estimator uses log returns)
def calc_R(X):
    return np.diff(np.log(X), axis=0)

########################################################################
#THESE FUNCTIONS ARE FOR THE DETECTOR
def G_hat(Omega):
    N=Omega.shape[0]
    P=np.eye(N)-(1/N)*np.ones((N, N))
    return P @ Omega @ P # Eqs. (13, 85)

def zero_sum_basis(N):
    Q = np.zeros((N, N - 1))
    for i in range(1, N):
        Q[:i, i - 1] = 1.0 / np.sqrt(i * (i + 1))
        Q[i, i - 1] = -i / np.sqrt(i * (i + 1))
    return Q

def null_basis(G, k):
    N = G.shape[0]
    Q = zero_sum_basis(N)
    lam, V = np.linalg.eigh(Q.T @ G @ Q)        # (N-1) eigenpairs, ascending
    J = (Q @ V[:, :k]).T                         # back to R^N: k x N, orthonormal, each sums to 0
    return J, lam[:k], lam

def alpha_hat(J, X, X_next, dt):
    r=np.log(X_next/X)
    return J@r/dt #Eq. (89): compute alpha A

def align(J_new, J_prev):
    if J_prev is None:
        return J_new
    C_tilde = J_new @ J_prev.T                   # Eq. (120): overlaps between new and old bases
    U, _, Vt = np.linalg.svd(C_tilde)
    C_hat = U @ Vt                               # nearest orthogonal matrix to C_tilde
    return C_hat.T @ J_new                       # rotate new basis onto the old one

def run_estimator(X, L=100, k=1, dt=1.0):
    """V&F detector"""
    T, N = X.shape[0] - 1, X.shape[1]
    logret = np.diff(np.log(X), axis=0)          # T x N returns (not log returns)
    A2 = np.full(T + 1, np.nan)
    ah = [None] * (T + 1)
    lams = [None] * (T + 1)
    J_values=[]
    J_prev, a_prev = None, None
    for t in range(L, T):                        # window = returns ending at t; predict step t -> t+1
        R = logret[t - L:t]                      # L x N
        Omega_hat = np.cov(R, rowvar=False, bias=True) / dt     # Eq. (93), 1/L normalisation
        G = G_hat(Omega_hat)                                     # Eq. (85)
        J, lam_k, lam_all = null_basis(G, k)                     # Eqs. (86)-(88)
        J = align(J, J_prev)                                     # Appendix A
        J_values.append(J)
        a = alpha_hat(J, X[t], X[t + 1], dt)                     # Eq. (89), uses J(t) and the NEXT step
        if a_prev is not None:
            A2[t + 1] = a @ a_prev                               # Eq. (90): sum_A alpha^A(t) alpha^A(t+1)
        ah[t + 1], lams[t + 1] = a, lam_all
        J_prev, a_prev = J, a
    return {"A2": A2, "alpha_hat": ah, "lambdas": lams, "J_hat": J_values}

#NOTE: When running on the Vf_series generated data, set use_cash to False!!!!!!!!!!!!!!!
def run_detector(R, L=250, refresh=5, n_dirs=1, use_cash=True, dt=1.0):
    """V&F detector on a (steps x n) array of returns.

    Parameters
    ----------
    R : ndarray, shape (T, n_in)
        Per-step returns (log returns per DEFINITIONS §1). One column per asset.
        If use_cash is True, a zero column is prepended, so n = n_in + 1;
        otherwise n = n_in and R must already contain cash if you want it.
    L : int
        Covariance window length. Window for step t is R[t-L+1 : t+1].
    refresh : int
        Re-estimate the calm portfolio(s) every `refresh` steps; reuse in between.
        refresh=1 re-estimates every step.
    n_dirs : int
        Number of calm portfolios k, 1 <= k <= n-1.
    use_cash : bool
        Prepend a column of zeros (cash, return 0) to R.
    dt : float
        Step size; divides Omega_hat and alpha_hat. Default 1 (per-step units).

    Returns
    -------
    dict with
    products : ndarray, shape (T,)
        A2[t+1] = sum_A alpha_A(t) * alpha_A(t+1)  (Eq. 90).
        nan for t <= L (no value until two consecutive alphas exist).
    port_ret : ndarray, shape (T,) if n_dirs == 1, else (T, n_dirs)
        alpha_A(t+1) = J^A(t) . r_{t+1} / dt  (Eq. 89). nan for t < L.
    weights : ndarray, shape (n_refresh, n) if n_dirs == 1, else (n_refresh, n_dirs, n)
        The portfolio(s) chosen at each refresh. Rows sum to 0, unit length.
        n_refresh = ceil((T - L) / refresh).
    refresh_times : ndarray, shape (n_refresh,)
        The step index t at which each row of `weights` was estimated.
    lambdas : ndarray, shape (n_refresh, n - 1)
        All eigenvalues of G_hat(t) on the zero-sum plane, ascending, at each refresh.
        Use to check the gap after the first n_dirs entries.
    """
    R = np.asarray(R, dtype=float)
    if use_cash:
        R = np.column_stack([np.zeros(len(R)), R])      # cash = column 0, return 0
    T, n = R.shape
    k = n_dirs
    assert 1 <= k <= n - 1 and T >= L + 2

    products = np.full(T, np.nan)
    port_ret = np.full((T, k), np.nan)
    weights, refresh_times = [], []

    J = J_prev = a_prev = None
    for t in range(L - 1, T - 1):
        if J is None or (t - (L - 1)) % refresh == 0:          # re-estimate every `refresh` steps
            W = R[t - L + 1: t + 1]
            Omega_hat = np.cov(W, rowvar=False, bias=True) / dt
            J, _, _ = null_basis(G_hat(Omega_hat), k)
            J = align(J, J_prev)
            J_prev = J
            weights.append(J.copy()); refresh_times.append(t)
        a = J @ R[t + 1] / dt                                   # Eq. 89
        if a_prev is not None:
            products[t + 1] = a @ a_prev                        # Eq. 90
        port_ret[t + 1] = a
        a_prev = a

    weights = np.array(weights)                                 # (n_refresh, k, n)
    if k == 1:
        weights = weights[:, 0, :]                              # (n_refresh, n) as the spec says
    return {"products": products, "port_ret": port_ret.squeeze(),
            "weights": weights, "refresh_times": np.array(refresh_times)}