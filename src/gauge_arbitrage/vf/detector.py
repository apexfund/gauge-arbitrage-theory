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
    r=(X_next-X)/X
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

def run_detector(R, L=100, k=1, dt=1.0):
    """V&F detector on a (steps x n) array of returns, any n >= 2.

    Returns dict:
        products : (T,)        A2[t+1] = sum_A alpha_A(t) * alpha_A(t+1)      (Eq. 90)
        alpha    : (T, k)      alpha_A(t+1) = J^A(t) . r_{t+1} / dt            (Eq. 89)
        weights  : (T, k, n)   the portfolio(s) J^A(t) chosen at each step
        lambdas  : (T, n-1)    all eigenvalues of G_hat(t), ascending
    Rows before the first defined output are nan.
    Window for Omega_hat(t) = R[t-L+1 : t+1]  (L returns ending at t inclusive).
    Include cash as a column of zeros if you want it in the universe.
    """
    R = np.asarray(R, dtype=float)
    if R.ndim != 2:
        raise ValueError("R must be a 2-D (steps x n) array of returns")
    T, n = R.shape
    if not 1 <= k <= n - 1:
        raise ValueError(f"k must be between 1 and n-1 = {n-1}")
    if T < L + 2:
        raise ValueError(f"need at least L+2 = {L+2} steps, got {T}")

    products = np.full(T, np.nan)
    alpha = np.full((T, k), np.nan)
    weights = np.full((T, k, n), np.nan)
    lambdas = np.full((T, n - 1), np.nan)

    J_prev, a_prev = None, None
    for t in range(L - 1, T - 1):                        # window ends at t; predict step t+1
        W = R[t - L + 1: t + 1]                          # (L, n)
        Omega_hat = np.cov(W, rowvar=False, bias=True) / dt   # Eq. (93)
        G = G_hat(Omega_hat)                             # Eqs. (13, 85)
        J, lam_k, lam_all = null_basis(G, k)             # Eqs. (86)-(88)
        J = align(J, J_prev)                             # Appendix A
        a = J @ R[t + 1] / dt                            # Eq. (89) on returns
        #   (identical to alpha_hat(J, X[t], X[t+1], dt) when R[t+1] = (X[t+1]-X[t])/X[t])
        if a_prev is not None:
            products[t + 1] = a @ a_prev                 # Eq. (90)
        alpha[t + 1], weights[t + 1], lambdas[t + 1] = a, J, lam_all
        J_prev, a_prev = J, a

    return {"products": products, "alpha": alpha, "weights": weights, "lambdas": lambdas}