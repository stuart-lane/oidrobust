import numpy as np

def newey_west(u_hat, Z, L, n):
    """
    Compute Newey-West variance estimator.
    
    Parameters
    ----------
    u_hat : array_like
        Residuals
    Z : array_like
        Instrument matrix
    L : int
        Number of lags
    n : int
        Sample size
        
    Returns
    -------
    array_like
        Variance matrix
    """
    V = Z.T @ np.diag(u_hat.flatten() ** 2) @ Z
    if L > 0:
        for l in range(1, L + 1):
            w_l = 1 - l / (L + 1)
            for t in range(l, n):
                gamma = w_l * u_hat[t] * u_hat[t - l]
                V += gamma * (Z[t:t+1].T @ Z[t-l:t-l+1] + 
                            Z[t-l:t-l+1].T @ Z[t:t+1])
    return V