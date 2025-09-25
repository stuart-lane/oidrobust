#' Score test for overidentification testing in linear IV models
#' 
#' @description
#' Implements various overidentification tests including the Sargan test,
#' Hansen J-test, and Kleibergen-Paap test for linear IV models. This function
#' implements the score test discussed in Windmeijer (2021) and Lane and Windmeijer (2025).
#' 
#' @param formula An object of class \code{formula} describing the IV regression model.
#'   Should include two parts separated by \code{|}, where the left part specifies
#'   the main equation and the right part specifies the instruments.
#' @param data A data frame containing the variables in the model.
#' @param y Numeric vector of the outcome variable.
#' @param X Matrix of endogenous regressors.
#' @param Z Matrix of instrumental variables.
#' @param W Matrix of included exogenous regressors. If NULL, an intercept is included.
#' @param weights Optional vector of weights.
#' @param method Character string specifying estimation method:
#'   \itemize{
#'     \item "2sls": Two-stage least squares (gives J-test or Sargan test)
#'     \item "liml": Limited information maximum likelihood (gives KP-test)
#'   }
#' @param errors Character string specifying variance estimator type:
#'   \itemize{
#'     \item "hom": Homoskedastic errors (gives Sargan test)
#'     \item "het": Heteroskedasticity-robust (White)
#'     \item "hac": HAC-robust (Newey-West)
#'     \item "cluster": Cluster-robust
#'   }
#' @param basmann Logical. If TRUE, uses Basmann's form of the test.
#' @param lags Specifies lags for Newey-West HAC estimator:
#'   \itemize{
#'     \item "rot": Rule of thumb (4*(n/100)^(2/9))
#'     \item "plug-in": Plug-in method ((4/3)*(n/100)^(1/4))
#'     \item numeric: User-specified number of lags
#'   }
#' @param cluster_var Vector specifying clusters for cluster-robust variance estimation.
#' @param small_sample Logical. If TRUE, applies small-sample correction.
#' @param no_constant Logical. If TRUE, suppresses the constant term.
#' 
#' @return An object of class "htest" containing:
#' \itemize{
#'   \item statistic: The score test statistic
#'   \item p.value: P-value based on chi-square distribution
#'   \item parameter: Degrees of freedom
#'   \item method: Description of test method used
#'   \item data.name: Description of input data
#'   \item coefficients: Estimated parameters
#' }
#' 
#' @details
#' The function implements overidentification tests that allow users to assess 
#' instrument validity when the number of instruments exceeds the number of 
#' endogenous regressors. It supports both traditional (Sargan) and robust 
#' (Hansen J, Kleibergen-Paap) test statistics.
#' 
#' The test statistic follows a chi-square distribution under the null hypothesis
#' of valid instruments, with degrees of freedom equal to the number of 
#' overidentifying restrictions (number of instruments minus number of 
#' endogenous regressors).
#' 
#' @examples
#' \dontrun{
#' # Using formula interface
#' data(mroz, package = "sampleSelection")
#' score_test(wage ~ educ | exper + city, data = mroz)
#' 
#' # Using matrix interface with robust errors
#' score_test(y = y, X = X, Z = Z, 
#'           method = "liml",
#'           errors = "hac",
#'           lags = "rot")
#'           
#' # Using clustered standard errors
#' score_test(wage ~ educ | exper + city, 
#'           data = mroz,
#'           errors = "cluster",
#'           cluster_var = "id")
#' }
#' 
#' @references
#' \itemize{
#'   \item Hansen, L. P. (1982). Large sample properties of generalized method of moments estimators. Econometrica, 1029-1054.
#'   \item Kleibergen, F., & Paap, R. (2006). Generalized reduced rank tests using the singular value decomposition. Journal of Econometrics, 133(1), 97-126.
#'   \item Lane, S., & Windmeijer, F. (2025). Overidentification testing with weak instruments and heteroskedasticity. Working paper
#'   \item Newey, W. K., & West, K. D. (1987). A simple, positive semi-definite, heteroskedasticity and autocorrelation consistent covariance matrix. Econometrica, 55(3), 703–708.
#'   \item Sargan, J. D. (1958). The estimation of economic relationships using instrumental variables. Econometrica, 393-415.
#'   \item Windmeijer, F. (2021). Testing underidentification in linear models, with applications to dynamic panel and asset pricing models. Journal of Econometrics, 105104.
#' }
#' 
#' @importFrom Formula as.Formula
#' @importFrom expm sqrtm
#' @importFrom stats model.frame model.matrix model.response terms pchisq
#' 
#' @export
### ============================================================================
### SCORE TEST FUNCTION
### ============================================================================

score_test <- function(formula = NULL, data = NULL, y = NULL, X = NULL, Z = NULL, W = NULL, 
                       weights = NULL, method = "2sls", errors = "hom", basmann = FALSE,
                       lags = "rot", cluster_var = NULL, small_sample = FALSE, no_constant = NULL) {
  
  if (!is.null(formula)) {
    if (!inherits(formula, "formula")) {
      stop("formula must be a valid R formula object")
    }
    if (is.null(data) || !is.data.frame(data)) {
      stop("data must be a provided data frame")
    }
    
    formula <- Formula::as.Formula(formula)
    if (length(formula)[2] != 2) {
      stop("Formula must have two right-hand sides")
    }
    
    mf <- model.frame(formula, data)
    y <- model.response(mf)
    
    X <- model.matrix(formula, data = mf, rhs = 1, intercept = TRUE)  
    Z <- model.matrix(formula, data = mf, rhs = 2, intercept = TRUE)  
    
    constant <- X[, "(Intercept)", drop = FALSE]
    
    X <- X[, -1, drop = FALSE]
    Z <- Z[, -1, drop = FALSE]
    
    exog_vars <- intersect(colnames(X), colnames(Z))
    
    W <- cbind(constant, X[, exog_vars, drop = FALSE])
    
    X <- X[, setdiff(colnames(X), exog_vars), drop = FALSE]
    Z <- Z[, setdiff(colnames(Z), exog_vars), drop = FALSE]
    
    if (!is.null(cluster_var)) {
      if (!cluster_var %in% colnames(data)) stop("Cluster variable not found in data")
      cluster_var <- data[[cluster_var]]
    }
  } else {
    
    if (is.null(y) || is.null(X) || is.null(Z)) {
      stop("When not using formula, y, X, and Z must be provided")
    }
    
    y <- as.matrix(y)
    X <- as.matrix(X)
    Z <- as.matrix(Z)
    
    if (is.null(W)) {
      
      W <- matrix(1, nrow = nrow(X), ncol = 1)
      
    } else {
      
      W <- as.matrix(W)
      
    }
    
    if (!is.null(cluster_var) && length(cluster_var) != nrow(X)) {
      stop("Cluster variable must have the same number of rows as the data")
    }
    
  }
  
  Z <- Z[, !duplicated(t(Z)), drop = FALSE]
  
  if (!is.null(weights)) {
    w_sqrt <- sqrt(weights)
    y <- y * w_sqrt
    X <- X * w_sqrt
    Z <- Z * w_sqrt
    W <- W * w_sqrt
  }
  
  if (!is.null(W)) {
    QR_W <- qr(W)
    Q_W <- qr.Q(QR_W)
    MW <- diag(nrow(W)) - Q_W %*% t(Q_W)
    y <- MW %*% y
    X <- MW %*% X
    Z <- MW %*% Z 
  }
  
  kx <- ncol(X)
  kz <- ncol(Z)
  if (kz - kx < 1) {
    stop("Model must be overidentified to conduct test")
  }
  
  n <- length(y)
  dof <- kz - kx
  
  pi_hat <- solve(t(Z) %*% Z) %*% t(Z) %*% X
  X_hat <- Z %*% pi_hat
  
  if (method == "2sls") {
    
    beta_hat <- solve(t(X_hat) %*% X) %*% t(X_hat) %*% y
    u_hat <- y - X %*% beta_hat
    M_X_hat <- diag(n) - X_hat %*% solve(t(X_hat) %*% X_hat) %*% t(X_hat)
    
  } else if (method == "liml") {
    
    Y <- cbind(y, X)
    sigY <- t(Y) %*% Y / n
    sigY12 <- sqrtm(sigY)
    YPZY <- t(Y) %*% Z %*% solve(t(Z) %*% Z) %*% t(Z) %*% Y
    eig_mat <- solve(sigY12) %*% YPZY %*% solve(sigY12) / n
    min_eig <- min(eigen(eig_mat)$values)
    
    beta_hat <- solve(t(X_hat) %*% X - min_eig * (t(X) %*% X)) %*% (t(X_hat) %*% y - min_eig * (t(X) %*% y))
    u_hat <- y - X %*% beta_hat
    M_u_hat <- diag(n) - u_hat %*% solve(t(u_hat) %*% u_hat) %*% t(u_hat)
    pi_liml <- solve(t(Z) %*% M_u_hat %*% Z) %*% t(Z) %*% M_u_hat %*% X
    X_hat_liml <- Z %*% pi_liml
    M_X_hat <- diag(n) - X_hat_liml %*% solve(t(X_hat_liml) %*% X_hat_liml) %*% t(X_hat_liml)
    
  } else {
    
    stop("Method must be either '2sls' or 'liml'")
    
  }
  
  Z2 <- Z[, (ncol(Z)-dof+1):ncol(Z)]
  
  MXZ2 <- M_X_hat %*% Z2
  
  if (!basmann) {
    residuals <- u_hat
  } else {
    MZ = diag(n) - (Z %*% solve(t(Z) %*% Z) %*% t(Z))
    residuals <- (MZ %*% y) - (MZ %*% X) %*% beta_hat
  }
  
  if (errors == "hom") {
    
    sigma2_hat <- (t(residuals) %*% residuals)
    sigma2_hat <- as.numeric(sigma2_hat)
    variance <- (t(MXZ2) %*% MXZ2) * sigma2_hat / n
    
  } else if (errors == "het") {
    
    residuals2 <- residuals^2
    variance <- t(MXZ2) %*% diag(as.vector(residuals2)) %*% MXZ2
    
  } else if (errors == "hac") {
    
    newey_west <- function(u_hat, Z, lags) {
      V <- t(Z) %*% diag(as.vector(u_hat^2)) %*% Z
      
      if (lags >0) {
        for (l in 1:lags) {
          w_l <- 1 - l / (L + 1)
          for (t in (l + 1):n) {
            V <- V + w_l * u_hat[t] * u_hat[t - l] * 
              (Z[t, ] %*% t(Z[t - l, ]) + Z[t - l, ] %*% t(Z[t, ]))
          }
        }
      } 
      
      return(V)
    }
    
    if (is.null(lags)) {
      stop("Select lags for Newey-West variance estimator")
    }
    if (lags == "rot") {
      L = floor(4 * (n / 100) ^ (2/9))
      print(paste0(L, " lags used for Newey-West estimator"))
    } else if (lags == "plug-in") {
      L = floor((4 / 3) * (n / 100) ^ (1/4))
      print(paste0(L, " lags used for Newey-West estimator"))
    } else {
      L = lags
    }
    
    variance <- newey_west(residuals, MXZ2, L)
    
  } else if (errors == "cluster") {
    
    if (is.null(cluster_var)) {
      stop("Select a cluster variable")
    }
    
    clusters <- unique(cluster_var)
    n_clusters <- length(clusters)
    
    variance <- matrix(0, ncol = dof, nrow = dof)
    
    for (cluster in clusters) {
      
      cluster_idx <- which(cluster_var == cluster)
      MXZ2_c <- as.matrix(MXZ2[cluster_idx, ])
      residuals_c <- as.matrix(u_hat[cluster_idx])
      variance <- variance + t(MXZ2_c) %*% residuals_c %*% t(residuals_c) %*% MXZ2_c
      
    }
    
    variance <- variance * (n_clusters / (n_clusters - 1))
    
  }
  
  score_stat <- t(u_hat) %*% MXZ2 %*% solve(variance) %*% t(MXZ2) %*% u_hat
  
  if (small_sample) {
    score_stat <- n * score_stat / (n - kz) 
  }
  
  p_value <- 1 - pchisq(score_stat, dof)
  
  test_statistic <- 
    if (method == "2sls") {
      if (errors == "het") {
        "J test for instrument validity (White heteroskedasticity-robust variance estimator)"
      } else if (errors == "hac") {
        paste0("J test for instrument validity (Newey-West HAC-robust variance estimator, ", L, " lags)")
      } else if (errors == "cluster") {
        paste0("J test for instrument validity (Cluster-robust variance estimator, ", n_clusters, " clusters)")
      } else {
        "2SLS-Sargan test instrument validity (Homoskedastic variance estimator)"
      }
    } else {
      if (errors == "het") {
        "KP test for instrument validity (White heteroskedasticity-robust variance estimator)"
      } else if (errors == "hac") {
        paste0("KP test for instrument validity (Newey-West HAC-robust variance estimator, ", L, " lags)")
      } else if (errors == "cluster") {
        paste0("KP test for instrument validity (Cluster-robust variance estimator, ", n_clusters, " clusters)")
      } else {
        "LIML-Sargan test instrument validity (Homoskedastic variance estimator)"
      }
    } 
  
  result <- list(
    statistic = c(score = as.numeric(score_stat)),
    p.value = p_value,
    parameter = c(dof = dof),
    method = paste(test_statistic),
    data.name = if (!is.null(formula)) deparse(substitute(data)) else "user-provided matrices",
    coefficients = beta_hat
  )
  class(result) <- "htest"
  
  return(result)
  
}
