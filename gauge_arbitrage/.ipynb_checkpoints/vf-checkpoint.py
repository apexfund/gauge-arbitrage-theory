import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.linalg import null_space

from .time_series import make_series, W_t

#\hat{\sigma}_mu^a calculation (right side of equation 108)
def calc_sigma_hat(sigma):
    #need to check this
    #sigma.shape[0] is number of rows (N)
    #sum over the rows (this makes the correct shape)
    return sigma-(1/(sigma.shape[0]))*np.sum(sigma, axis=0)

#write these functions tomorrow (discuss with apex pod)
def calc_J(sigma):
    #J dot 1 = 0 (sum of components is 0), J dot sigma^a = 0 for all a (a is column index)
    #>=N-1-d (more if there are linearly dependent columns in sigma)
    return null_space(np.column_stack([sigma, np.ones(sigma.shape[0])]).T)

#FOR V&F SECTION 4.3 ONLY
#equation 110
def calc_alpha_43(sigma, beta, alpha_A):
    return (1/sigma.shape[0])*np.sum(np.matmul(beta, calc_sigma_hat(sigma).T))-np.sum(np.matmul(alpha_A, calc_J(sigma).T))

#FOR V&F SECTION 4.3 ONLY
#equation 108 (left)
def calc_alpha_mu_43(mu, sigma, beta, alpha_A):
    if mu==0:
        return 0
    else:
        return calc_alpha_43(sigma, beta, alpha_A)+np.sum(beta @ calc_sigma_hat(sigma)[mu])+np.sum(alpha_A @ calc_J(sigma)[mu])

#v&f equation 2 SDE solution (equation 106)
#alpha_mu is a 1d array, sigma (sigma_mu^a) is 2d array (mu: 1 to N=row, a: 1 to d=column)
#pls tell me if theres a way to make this more efficient :skull:
@make_series
def X_mu(T, n, alpha_mu, sigma, X_0, mu, t=None):
    B=np.array([W_t(T, n) for i in range(sigma[mu].size)]) #define array of B (B(t)) in the formula.
    #this should create d=sigma[mu].size (# cols) Brownian motions, each of size n+1
    #We broadcast, NOT matrix multiply (that's what sigma[mu]*B.T does) and sum along d
    return X_0*np.exp((alpha_mu[mu]-0.5*np.sum(np.square(sigma[mu])))*t + np.sum(sigma[mu]*B.T, axis=1))

#FOR V&F SECTION 4.3 ONLY (same as the last function but i plugged the other section 4.3 functions in)
@make_series
def X_mu_43(T, n, sigma, beta, alpha_A, X_0, mu, t=None):
    B=np.array([W_t(T, n) for i in range(sigma[mu].size)]) 
    return X_0*np.exp((calc_alpha_mu_43(mu, sigma, beta, alpha_A)-
                       0.5*np.sum(np.square(sigma[mu])))*t + np.sum(sigma[mu]*B.T, axis=1))

if __name__ == "__main__":
    pass