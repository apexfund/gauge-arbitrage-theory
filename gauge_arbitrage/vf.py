import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.linalg import null_space

from .time_series import Series

########################################################################
#THESE FUNCTIONS ARE FOR MAKING THE RANDOM WALKS
def calc_sigma_hat(sigma):
    """\hat{\sigma}_mu^a calculation (right side of equation 108)"""
    #need to check this
    #sigma.shape[0] is number of rows (N)
    #sum over the rows (this makes the correct shape)
    return sigma-(1/(sigma.shape[0]))*np.sum(sigma, axis=0)

def calc_J(sigma):
    """Calculates J based on equation 5"""
    #J dot 1 = 0 (sum of components is 0), J dot sigma^a = 0 for all a (a is column index)
    #>=N-1-d (more if there are linearly dependent columns in sigma)
    return null_space(np.column_stack([sigma, np.ones(sigma.shape[0])]).T)

#FOR V&F SECTION 4.3 ONLY
#equation 110
def calc_alpha_43(sigma, beta, alpha_A):
    """Calculates alpha for the section 4.3 simulation"""
    return (1/sigma.shape[0])*np.sum(np.matmul(beta, calc_sigma_hat(sigma).T))-np.sum(np.matmul(alpha_A, calc_J(sigma).T))

#FOR V&F SECTION 4.3 ONLY
#equation 108 (left)
def calc_alpha_mu_43(sigma, beta, alpha_A):
    """Calculates alpha_mu for the section 4.3 simulation"""
    alpha_mu=calc_alpha_43(sigma, beta, alpha_A)+(beta @ calc_sigma_hat(sigma).T)+(alpha_A @ calc_J(sigma).T)
    alpha_mu-=alpha_mu[0]
    return alpha_mu

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
    logret = np.diff(np.log(X), axis=0)          # T x N log-returns
    A2 = np.full(T + 1, np.nan)
    ah = [None] * (T + 1)
    lams = [None] * (T + 1)
    J_prev, a_prev = None, None
    for t in range(L, T):                        # window = returns ending at t; predict step t -> t+1
        R = logret[t - L:t]                      # L x N
        Omega_hat = np.cov(R, rowvar=False, bias=True) / dt     # Eq. (93), 1/L normalisation
        G = G_hat(Omega_hat)                                     # Eq. (85)
        J, lam_k, lam_all = null_basis(G, k)                     # Eqs. (86)-(88)
        J = align(J, J_prev)                                     # Appendix A
        a = alpha_hat(J, X[t], X[t + 1], dt)                     # Eq. (89), uses J(t) and the NEXT step
        if a_prev is not None:
            A2[t + 1] = a @ a_prev                               # Eq. (90): sum_A alpha^A(t) alpha^A(t+1)
        ah[t + 1], lams[t + 1] = a, lam_all
        J_prev, a_prev = J, a
    return {"A2": A2, "alpha_hat": ah, "lambdas": lams}

########################################################################
#Subclass for the series in V&F (Eq. (106))
class vf_series(Series):
    #alpha_mu is a 1d array, sigma (sigma_mu^a) is 2d array (mu: 1 to N=row, a: 1 to d=column)
    #pls tell me if theres a way to make this more efficient :skull:
    def X_mu(self, alpha_mu, sigma, X_0, mu):
        """Solution to equation 106"""
        B=np.array([self.bm() for i in range(sigma[mu].size)]) #define array of B (B(t)) in the formula.
        #this should create d=sigma[mu].size (# cols) Brownian motions, each of size n+1
        #We broadcast, NOT matrix multiply (that's what sigma[mu]*B.T does) and sum along d
        return X_0*np.exp((alpha_mu[mu]-0.5*np.sum(np.square(sigma[mu])))*self.t + np.sum(sigma[mu]*B.T, axis=1))

    #FOR V&F SECTION 4.3 ONLY (same as the last function but i plugged the other section 4.3 functions in)
    def X_mu_43(self, sigma, beta, alpha_A, X_0, Z):
        """Solution to equation 106 for section 4.3 simulation"""
        B=np.array([self.vf_bm(Z, i) for i in range(sigma.shape[1])]) 
        return X_0*np.exp((calc_alpha_mu_43(sigma, beta, alpha_A)-
                           0.5*np.sum(np.square(sigma), axis=1)) * self.t[:, None] + (sigma @ B).T)
    
    def vf_bm(self, Z, i):
        dW=np.sqrt(self.dt) * Z[:, i]
        W = np.concatenate([[0], np.cumsum(dW)])
        return W

if __name__ == "__main__":
    pass