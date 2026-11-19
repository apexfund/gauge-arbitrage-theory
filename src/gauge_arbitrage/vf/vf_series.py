import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.linalg import null_space

from gauge_arbitrage.time_series import Series

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
#Subclass for the series in V&F (Eq. (106))
class Vf_series(Series):
    #alpha_mu is a 1d array, sigma (sigma_mu^a) is 2d array (mu: 1 to N=row, a: 1 to d=column)
    #pls tell me if theres a way to make this more efficient :skull:
    def X_mu(self, alpha_mu, sigma, X_0, mu):
        """Solution to equation 106"""
        B=np.array([self.bm() for i in range(sigma[mu].size)]) #define array of B (B(t)) in the formula.
        #this should create d=sigma[mu].size (# cols) Brownian motions, each of size n+1
        #We broadcast, NOT matrix multiply (that's what sigma[mu]*B.T does) and sum along d
        return X_0*np.exp((alpha_mu[mu]-0.5*np.sum(np.square(sigma[mu])))*self.times + np.sum(sigma[mu]*B.T, axis=1))

    #FOR V&F SECTION 4.3 ONLY (same as the last function but i plugged the other section 4.3 functions in)
    def X_mu_43(self, sigma, beta, alpha_A, X_0, Z):
        """Solution to equation 106 for section 4.3 simulation"""
        B=np.array([self.vf_bm(Z, i) for i in range(sigma.shape[1])]) 
        return X_0*np.exp((calc_alpha_mu_43(sigma, beta, alpha_A)-
                           0.5*np.sum(np.square(sigma), axis=1)) * self.times[:, None] + (sigma @ B).T)
    
    def vf_bm(self, Z, i):
        dW=np.sqrt(self.dt) * Z[:, i]
        W = np.concatenate([[0], np.cumsum(dW)])
        return W
