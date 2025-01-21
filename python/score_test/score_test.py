"""
Package for overidentification testing in linear IV models
"""

import numpy as np
from scipy import stats, linalg
from scipy.linalg import sqrtm
import pandas as pd
from patsy import dmatrices
from typing import Optional, Union, Dict, Any
from dataclasses import dataclass

@dataclass
class TestResults:
    """Container for test results"""
    statistic: float
    p_value: float
    df: int
    coefficients: np.ndarray
    method: str

    def __str__(self) -> str:
        return (
            f"Overidentification test results\n"
            f"-------------------------------\n"
            f"Method: {self.method}\n"
            f"Statistic: {self.statistic}\n"
            f"p-value: {self.p_value}\n"
            f"Degrees of freedom: {self.df}\n"
        )
    
class VarianceEstimator:
    """Standard variance covariance matrix estimator"""

    def compute(self, residuals: np.ndarray, MXZ2: np.ndarray, **kwargs) -> np.ndarray:
        raise NotImplementedError
    
class HomoskedasticVariance(VarianceEstimator):
    """Homoskedastic variance estimator"""

    def compute(self, residuals: np.ndarray, MXZ2: np.ndarray, **kwargs) -> np.ndarray:
        n = len(residuals)
        sigma2_hat = np.sum(residuals ** 2) / n
        return (MXZ2.T @ MXZ2) * sigma2_hat
    
class HeteroskedasticVariance(VarianceEstimator):
    """White heteroskedasticity-robust variance estimator"""

    def compute(self, residuals: np.ndarray, MXZ2: np.ndarray, **kwargs) -> np.ndarray:
        residuals2 = residuals ** 2
        return MXZ2.T @ np.diag(residuals2.flatten()) @ MXZ2
    
class NeweyWestVarianceEstimator(VarianceEstimator):
    """Newey-West Heteroskedasticity and autocorrelation robust variance estimator"""

    def _get_lag_length(self, n: int, method: str) -> int:
        if method == "rot":
            return int(4 * (n / 100) ** (2 / 9))
        elif method == "plug-in":
            return int((4 / 3) * (n / 100) ** (1 / 4))
        else:
            try:
                return int(method)
            except ValueError:
                raise ValueError("Please enter 'rot', 'plug-in' or a positive integer for the `lags` argument")
            
    def compute(self, residuals: np.ndarray, MXZ2: np.ndarray,
                lags: Union[str, int] = "rot", **kwargs) -> np.ndarray:
        n = len(residuals)
        V = MXZ2.T @ np.diag(residuals.flatten() ** 2) @ MXZ2
        L = self._get_lag_length(n, lags)
        
        if L > 0:
            for l in range(1, L + 1):
                w_l = 1 - l / (L + 1)
                for t in range(l, n):
                    gamma = w_l * residuals[t] * residuals[t - l]
                    V += gamma * (MXZ2[t:t+1].T @ MXZ2[t-l:t-l+1] + 
                                MXZ2[t-l:t-l+1].T @ MXZ2[t:t+1])
        return V
    
class ClusterVariance(VarianceEstimator):
    """Cluster robust variance estimator"""

    def compute(self, residuals: np.ndarray, MXZ2: np.ndarray,
                cluster_var: np.ndarray, **kwargs) -> np.ndarray:
        if cluster_var is None:
            raise ValueError("Cluster variable must be provided for cluster-robust variance errors")
        
        n = len(residuals)
        clusters = np.unique(cluster_var)
        n_clusters = len(clusters)
        df = MXZ2.shape[1]
        V = np.zeros((df, df))

        for cluster in clusters:
            cluster_idx = np.where(cluster_var == cluster)[0]
            MXZ2_c = MXZ2[cluster_idx]
            residuals_c = residuals[cluster_idx]
            V += MXZ2_c.T @ residuals_c @ residuals_c.T @ MXZ2_c

        V += V * n_clusters / (n_clusters - 1)
        return V

class ScoreTest:
    """Class for conducting overidentifcation tests in linear IV models"""

    def __init__(self):
        self.variance_estimators = {
            "hom": HomoskedasticVariance(),
            "het": HeteroskedasticVariance(),
            "hac": NeweyWestVarianceEstimator(),
            "cluster": ClusterVariance()
        }

    def _prepare_data(self,
                    y: Optional[np.ndarray] = None,
                    X: Optional[np.ndarray] = None,
                    Z: Optional[np.ndarray] = None,
                    W: Optional[np.ndarray] = None,
                    weights: Optional[np.ndarray] = None,
                    no_constant: bool = False) -> tuple:
        """
        Prepare data for overidentification testing.
        """

        # Basic error check
        if y is None or X is None or Z is None:
            raise ValueError("y, X and Z must be provided.")
        
        # Ensure correct shapes
        y = np.asarray(y).reshape(-1, 1)
        X = np.asarray(X).reshape(-1, 1) if X.ndim == 1 else np.asarray(X)
        Z = np.asarray(Z)

        # Handle constant
        if not no_constant:
            if W is None:
                W = np.ones((len(y), 1))
            else:
                W = np.asarray(W)
                has_constant = np.any(np.all(np.abs(W - np.ones((W.shape[0], 1))) < 1e-8, axis=0))
                if not has_constant:
                    W = np.column_stack([np.ones(len(y)), W])

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

        return y, X, Z, W

    def _partial_out(self, y: np.ndarray, X: np.ndarray, Z: np.ndarray,
                    W = Optional[np.ndarray]) -> tuple:
        """Partial out any exogenous regressors, including constant"""
        if W is not None:
            Q_W = np.linalg.qr(W)[0]
            M_W = np.eye(len(y)) - Q_W @ Q_W.T
            y = M_W @ y
            X = M_W @ X 
            Z = M_W @ Z
        return y, X, Z
    
    def compute_2sls(self, y: np.ndarray, X: np.ndarray, Z: np.ndarray) -> tuple:
        """Compute the 2SLS estimator"""
        pi_hat = np.linalg.pinv(Z.T @ Z) @ (Z.T @ X)
        X_hat = Z @ pi_hat
        beta_hat = np.linalg.pinv(X_hat.T @ X_hat) @ (X_hat.T @ y)
        u_hat = y - X @ beta_hat
        M_X_hat = np.eye(len(y)) - X_hat @ np.linalg.pinv(X_hat.T @ X_hat) @ X_hat.T
        return beta_hat, u_hat, M_X_hat
        
    def compute_liml(self, y: np.ndarray, X: np.ndarray, Z: np.ndarray) -> tuple:
        """Compute the LIML estimator"""
        n = len(y)
        Y = np.column_stack([y, X])
        sigY = Y.T @ Y / n
        sigY12 = sqrtm(sigY)
        YPZY = Y.T @ Z @ np.linalg.pinv(Z.T @ Z) @ Z.T @ Y
        eig_mat = np.linalg.solve(sigY12, YPZY @ np.linalg.solve(sigY12.T, np.eye(sigY12.shape[0]))) / n
        min_eig = min(np.linalg.eigvals(eig_mat).real)
        
        X_hat = Z @ np.linalg.pinv(Z.T @ Z) @ Z.T @ X
        beta_hat = np.linalg.pinv(X_hat.T @ X - min_eig * (X.T @ X)) @ (X_hat.T @ y - min_eig * (X.T @ y))
        u_hat = y - X @ beta_hat
        
        M_u_hat = np.eye(n) - u_hat @ np.linalg.pinv(u_hat.T @ u_hat) @ u_hat.T
        pi_liml = np.linalg.pinv(Z.T @ M_u_hat @ Z) @ Z.T @ M_u_hat @ X
        X_hat_liml = Z @ pi_liml
        M_X_hat = np.eye(n) - X_hat_liml @ np.linalg.pinv(X_hat_liml.T @ X_hat_liml) @ X_hat_liml.T
        
        return beta_hat, u_hat, M_X_hat
        
    def score_test(self,
                    formula: Optional[str] = None,
                    data: Optional[pd.DataFrame] = None,
                    y: Optional[np.ndarray] = None,
                    X: Optional[np.ndarray] = None,
                    Z: Optional[np.ndarray] = None,
                    W: Optional[np.ndarray] = None,
                    weights: Optional[np.ndarray] = None,
                    method: str = "2sls",
                    errors: str = "hom",
                    basmann: bool = False,
                    lags: Union[str, int] = "rot",
                    cluster_var: Optional[np.ndarray] = None,
                    small_sample: bool = False,
                    no_constant: bool = False) -> TestResults:

        """Conduct overidentification test in linear IV models
        
        ---- Paramaters ----

        // Formula parameters
        formula : str, optional
            Patsy formula specification
        data : pd.DataFrame, optional
            Pandas dataframe for formula

        // User-defined matrix parameters
        y : array_like, optional
            Dependent variable
        X : array_like, optional
            Endogenous regressors
        Z : array_like, optional
            Excluded instruments
        W : arrray_like, optional,
            Included exogenous regressors
        
        // Other parameters
        weights : array_like, optional
            Weights for computing regression coefficients
        method : str, optional
            Method used to calculate regression coefficient ('2sls' or 'liml')
            Default is '2sls'
        errors : str, optional 
            Structure of the variance matrix ('hom', 'het', 'hac' or 'cluster')
            Default is 'hom'
        basmann : bool, optional
            Use Basmann version of the test statistic
            Default is False
        lags : str or int, optional
            Number of lags for HAC variance estimator ('rot', 'plug-in' or user-specified positive integer)
            Default is 'rot'
        cluser_var : array_like, optional
            Clustering variable for cluster variance estimator
        small_sample : bool, optional
            Apply small-sample variance estimator correction
            Default is False
        no_constant : bool, optional
            Do not add constant to regression
            Default is False            

        ---- Returns ---- 

        TestResults
            Test results container
        """

        y, X, Z, W = self._prepare_data(
            y = y,
            X = X,
            Z = Z,
            W = W,
            weights = weights,
            no_constant = no_constant
        )

        # Ensure overidentification for computability of test statistic
        n, kx = X.shape[0], X.shape[1]
        kz = Z.shape[1]
        if kz - kx < 1:
            raise ValueError("Model must be overidentified.")
        df = kz - kx

        # Partial out any included exogenous regressors
        y, X, Z = self._partial_out(y, X, Z, W)

        if method.lower() == "2sls":
            beta_hat, u_hat, M_X_hat = self.compute_2sls(y, X, Z)
        elif method.lower() == "liml":
            beta_hat, u_hat, M_X_hat = self.compute_liml(y, X, Z)
        else:
            raise ValueError("Method specified must be either '2sls' or 'liml'.")
        
        Z2 = Z[:, -df:]
        MXZ2 = M_X_hat @ Z2

        if not basmann:
            residuals = u_hat
        else:
            MZ = np.eye(n) - Z @ np.linalg.pin(Z.T @ Z) @ Z.T
            residuals = MZ @ y - MZ @ (X @ beta_hat)

        variance_matrix = self.variance_estimators.get(errors.lower())
        if variance_matrix is None:
            raise ValueError(f"Unknown error structure: {errors}. Please select from 'hom', 'het', 'hac' or 'cluster'.")
        
        variance = variance_matrix.compute(
            residuals = residuals,
            MXZ2 = MXZ2,
            lags = lags,
            cluster_var = cluster_var
        )

        # Compute the score statistic
        score_stat = float(residuals.T @ MXZ2 @ np.linalg.pinv(variance) @ MXZ2.T @ residuals)

        if small_sample:
            score_stat = n * score_stat / (n - kz)


        p_value = 1 - stats.chi2.cdf(score_stat, df)

        return TestResults(
            statistic=score_stat,
            p_value=p_value,
            df=df,
            coefficients=beta_hat,
            method=f"{'2SLS' if method.lower() == '2sls' else 'LIML'}"
        )