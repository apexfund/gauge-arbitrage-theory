import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from functools import wraps #for the decorator
import inspect #for the decorator
from scipy.linalg import null_space

class Series:  
    def __init__(self, T, n, seed=0):
        self.T=T
        self.n=n
        self.rng=np.random.default_rng(seed)
    
    @property
    def dt(self):
        return self.T/self.n
    
    @property
    def t(self):
        return np.linspace(0.0, self.T, self.n+1)
    
    def bm(self):
        """Brownian motion W(t)"""
        dW=np.sqrt(self.dt) * self.rng.standard_normal(self.n) # increments ~ N(0, dt) (required for Brownian motion)
        W = np.concatenate([[0], np.cumsum(dW)]) #running total and concatenate to 0
        return W
    
    def random_walk(self):
        """Random Walk"""
        return self.bm()

    def mr_noise(self,rho): 
        """create mean-reverting noise"""
        e = np.sqrt(self.dt) * rng.standard_normal(self.n)
        #create mean-reverting noise
        a=np.array([0])
        for i in range(n):
            a=np.append(a, rho*a[i]+e[i]) #def of mean reversion
        return a

    def abm(self, alpha, sigma, X_0):
        """Create arithmetic Brownian motion"""
        return X_0+(alpha*self.t)+(sigma*self.bm())

    def gbm(self, alpha, sigma, X_0):
        """Create geometric Brownian motion"""
        return X_0*np.exp((alpha-(0.5*sigma**2))*self.t+(sigma*self.bm()))


if __name__ == "__main__":
    pass