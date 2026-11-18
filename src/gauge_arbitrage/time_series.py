import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from functools import wraps #for the decorator
import inspect #for the decorator
from scipy.linalg import null_space

class Series:  
    """Creates time series for analysis using NumPy.
    Example:
    >>> from gauge_arbitrage import time_series

    >>> A=time_series.Series(100, 100, seed=0)
    >>> fig, ax=plt.subplots(1,1)

    >>> ax.plot(A.times, A.bm())              # Brownian Motion (random walk)
    >>> ax.plot(A.times, A.abm(0.5,0.5,0))    # Arithmetic Brownian Motion
    >>> ax.plot(A.times, A.gbm(0.08, 0.2, 1)) # Geometric Brownian Motion"""
    def __init__(self, horizon, T, seed=0):
        self.horizon=horizon
        self.T=T
        self.rng=np.random.default_rng(seed)
    
    @property
    def dt(self):
        return self.horizon/self.T
    
    @property
    def times(self):
        return np.linspace(0.0, self.horizon, self.T+1)
    
    def bm(self):
        """Brownian motion W(t)"""
        dW=np.sqrt(self.dt) * self.rng.standard_normal(self.T) # increments ~ N(0, dt) (required for Brownian motion)
        W = np.concatenate([[0], np.cumsum(dW)]) #running total and concatenate to 0
        return W
    
    def random_walk(self):
        """Random Walk"""
        return self.bm()

    def mr_noise(self,rho): 
        """create mean-reverting noise"""
        e = np.sqrt(self.dt) * rng.standard_normal(self.T)
        #create mean-reverting noise
        a=np.array([0])
        for i in range(self.T):
            a=np.append(a, rho*a[i]+e[i]) #def of mean reversion
        return a

    def abm(self, alpha, sigma, X_0):
        """Create arithmetic Brownian motion"""
        return X_0+(alpha*self.times)+(sigma*self.bm())

    def gbm(self, alpha, sigma, X_0):
        """Create geometric Brownian motion"""
        return X_0*np.exp((alpha-(0.5*sigma**2))*self.times+(sigma*self.bm()))


if __name__ == "__main__":
    pass