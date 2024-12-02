import numpy as np
from scipy import stats, linalg
from scipy.linalg import sqrtm
import pandas as pd
from patsy import dmatrices

def score_test(formula=None, data=None, y=None, X=None, Z=None, W=None, 
               weights=None, method="2sls", errors="hom", basmann=False, 
               lags="rot", cluster_var=None, small_sample=False, no_constant=False):
    """
    Score test for overidentification testing in linear IV models.
    """
    # Handle formula interface
    if formula is not None:
        if data is None:
            raise ValueError("Data must be provided when using formula")
        if not isinstance(data, pd.DataFrame):
            raise ValueError("data must be a DataFrame")
            
        # Parse formula and create matrices
        y, X = dmatrices(formula.split('|')[0], data)
        _, Z = dmatrices(formula.split('|')[1], data)
        
        if not no_constant:
            if W is None:
                W = np.ones((len(y), 1))
            elif not np.all(W[:, 0] == 1):
                W = np.column_stack([np.ones(len(y)), W])

    # Basic error checks
    if y is None or X is None or Z is None:
        raise ValueError("When not using formula, y, X, and Z must be provided")

    # Convert inputs to numpy arrays
    y = np.asarray(y).reshape(-1, 1)
    X = np.asarray(X).reshape(-1, 1)
    Z = np.asarray(Z)
    
    # Remove duplicate instruments
    _, unique_cols = np.unique(Z, axis=1, return_index=True)
    Z = Z[:, np.sort(unique_cols)]

    # Apply weights if provided
    if weights is not None:
        w_sqrt = np.sqrt(weights)
        y = y * w_sqrt.reshape(-1, 1)
        X = X * w_sqrt.reshape(-1, 1)
        Z = Z * w_sqrt.reshape(-1, 1)
        if W is not None:
            W = W * w_sqrt.reshape(-1, 1)

    # Partial out exogenous regressors
    if W is not None:
        Q_W = np.linalg.qr(W)[0]
        M_W = np.eye(len(W)) - Q_W @ Q_W.T
        y = M_W @ y
        X = M_W @ X
        Z = M_W @ Z

    # Check for overidentification
    n, kx = X.shape[0], X.shape[1]
    kz = Z.shape[1]
    if kz - kx < 1:
        raise ValueError("Model must be overidentified to conduct test")
    df = kz - kx

    # First stage
    pi_hat = np.linalg.pinv(Z.T @ Z) @ Z.T @ X
    X_hat = Z @ pi_hat

    if method == "2sls":
        # 2SLS estimation
        beta_hat = np.linalg.pinv(X_hat.T @ X) @ X_hat.T @ y
        u_hat = y - X @ beta_hat
        M_X_hat = np.eye(n) - X_hat @ np.linalg.pinv(X_hat.T @ X_hat) @ X_hat.T
    
    elif method == "liml":
        # LIML estimation
        Y = np.column_stack([y, X])
        sigY = Y.T @ Y / n
        sigY12 = sqrtm(sigY)
        YPZY = Y.T @ Z @ np.linalg.pinv(Z.T @ Z) @ Z.T @ Y
        eig_mat = np.linalg.solve(sigY12, YPZY @ np.linalg.solve(sigY12.T, np.eye(sigY12.shape[0]))) / n
        min_eig = min(np.linalg.eigvals(eig_mat))
        
        beta_hat = np.linalg.pinv(X_hat.T @ X - min_eig * (X.T @ X)) @ (X_hat.T @ y - min_eig * (X.T @ y))
        u_hat = y - X @ beta_hat
        M_u_hat = np.eye(n) - u_hat @ np.linalg.pinv(u_hat.T @ u_hat) @ u_hat.T
        pi_liml = np.linalg.pinv(Z.T @ M_u_hat @ Z) @ Z.T @ M_u_hat @ X
        X_hat_liml = Z @ pi_liml
        M_X_hat = np.eye(n) - X_hat_liml @ np.linalg.pinv(X_hat_liml.T @ X_hat_liml) @ X_hat_liml.T

    # Extract overidentifying instruments
    Z2 = Z[:, -df:]
    MXZ2 = M_X_hat @ Z2

    # Compute residuals based on test type
    if not basmann:
        residuals = u_hat
    else:
        MZ = np.eye(n) - Z @ np.linalg.pinv(Z.T @ Z) @ Z.T
        residuals = MZ @ y - MZ @ X @ beta_hat

    # Compute variance based on error structure
    if errors == "hom":
        sigma2_hat = residuals.T @ residuals
        variance = (MXZ2.T @ MXZ2) * sigma2_hat / n
    
    elif errors == "het":
        residuals2 = residuals ** 2
        variance = MXZ2.T @ np.diag(residuals2.flatten()) @ MXZ2
    
    elif errors == "hac":
        if isinstance(lags, str):
            if lags == "rot":
                L = int(4 * (n / 100) ** (2/9))
            elif lags == "plug-in":
                L = int((4/3) * (n / 100) ** (1/4))
        else:
            L = int(lags)
        
        variance = newey_west(residuals, MXZ2, L)
    
    elif errors == "cluster":
        if cluster_var is None:
            raise ValueError("Cluster variable must be provided for cluster-robust errors")
        
        clusters = np.unique(cluster_var)
        n_clusters = len(clusters)
        variance = np.zeros((df, df))
        
        for cluster in clusters:
            cluster_idx = np.where(cluster_var == cluster)[0]
            MXZ2_c = MXZ2[cluster_idx]
            residuals_c = residuals[cluster_idx]
            variance += MXZ2_c.T @ residuals_c @ residuals_c.T @ MXZ2_c
        
        variance *= n_clusters / (n_clusters - 1)

    # Compute score statistic
    score_stat = float(y.T @ MXZ2 @ np.linalg.pinv(variance) @ MXZ2.T @ y)
    
    if small_sample:
        score_stat = n * score_stat / (n - kz)
    
    p_value = 1 - stats.chi2.cdf(score_stat, df)
    
    return {
        'statistic': score_stat,
        'p_value': p_value,
        'df': df,
        'coefficients': beta_hat,
        'method': f"{'2SLS' if method == '2sls' else 'LIML'}-Score test ({errors} errors)"
    }