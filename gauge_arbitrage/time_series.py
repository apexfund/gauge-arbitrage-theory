import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from functools import wraps #for the decorator
import inspect #for the decorator
from scipy.linalg import null_space

rng = np.random.default_rng(seed=0)

def make_series(f):
    sig = inspect.signature(f) #basically this is just so I can use this on any function with parameters T and n
    wants_t = "t" in sig.parameters #in the case the function requires t (abm and gbm)
    @wraps(f)
    def wrapper(*args, **kwargs):
        bound = sig.bind(*args, **kwargs)  
        bound.apply_defaults()              
        #get the T and n
        T = bound.arguments["T"]
        n = bound.arguments["n"]
        t = np.linspace(0, T, n + 1) #create the array of times (intervals of dt)
        if wants_t and "t" not in kwargs:
            kwargs["t"] = t
        return pd.DataFrame({'Value':f(*args, **kwargs), 'Time':t}) #return dataframe
    return wrapper

#Brownian motion W(t)
def W_t(T, n):
    dt = T/n
    dW = np.sqrt(dt) * rng.standard_normal(n) # increments ~ N(0, dt) (required for Brownian motion)
    W = np.concatenate([[0], np.cumsum(dW)]) #running total and concatenate to 0
    return W

#Random walk
@make_series
def random_walk(T,n):
    return W_t(T, n)

#create mean-reverting noise
def mr_noise(T,n,rho):
    dt = T/n 
    t = np.linspace(0, T, n + 1)
    e = np.sqrt(dt) * rng.standard_normal(n)
    #create mean-reverting noise
    a=np.array([0])
    for i in range(n):
        a=np.append(a, rho*a[i]+e[i]) #def of mean reversion
    return a

#arithmetic brownian motion
@make_series
def abm(T, n, alpha, sigma, X_0, t=None):
    return X_0+(alpha*t)+(sigma*W_t(T, n))

#geometric brownian motion (log normal random walks)
@make_series
def gbm(T, n, alpha, sigma, X_0, t=None):
    return X_0*np.exp((alpha-(0.5*sigma**2))*t+(sigma*W_t(T,n)))

if __name__ == "__main__":
    pass