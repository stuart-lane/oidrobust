Function arguments and basic usage
================
Stuart Lane
2025-09-25

The `score_test()` function takes on the following arguments:

    score_test(formula = NULL, data = NULL, y = NULL, X = NULL, Z = NULL, W = NULL, 
               weights = NULL, method = "2sls", errors = "hom", basmann = FALSE,
               lags = "rot", cluster_var = NULL, small_sample = FALSE, no_constant = NULL)

### Detailed description of the arguments

#### Core arguments

- `y`:
  - Outcome of interest
  - Must be specified
  - Must be a vector
- `X`:
  - Endogenous regressor(s)
  - Must be specified
  - Must be a vector or matrix
  - Must have the same length as `y`
- `Z`:
  - Instruments
  - Must be specified
  - Must be a matrix
  - Must have the same length as `y`
  - Must have more columns than `X` so that the model is overidentified
- `W`:
  - Included exogenous variable(s)
  - Optional, default is `NULL`
  - If included, it must be a vector or matrix
  - If included, it must have the same length as `y`
- `weights`:
  - Specifies regression weights
  - Optional, default is `NULL`
  - If included, it must be a vector
  - If included, it must have the same length as `y`
- `method`:
  - Specifies using 2SLS or LIML for the score test
  - Default is `"2sls"`
  - Must be either `"2sls"` or `"liml"`

#### Variance arguments

- `errors`:
  - Specifies the type of covariance-variance matrix estimator used
  - Optional, efault is `"hom"` (homoskedastic)
  - If specified, must be one of these options:
    - `"hom"` (homoskedastic variance estimator)
    - `"het"` (White heteroskedasticity-robust variance estimator)
    - `"hac"` (Newey-West heteroskedasticity- and autocorrelation-robust
      variance estimator)
    - `"cluster"` (cluster-robust variance estimator)

If `errors = "hac"`, then there is an additional optional argument:

- `lags`:
  - Specifies the number of lags for the Newey-West variance estimator
  - Optional, default is `"rot"` (rule-of-thumb)
  - Available preset options to automatically calculate number of lags
    are:
    - `"rot"` (rule-of-thumb)
    - `"plug-in"` (plug-in)
  - Alternatively, the user can specify any non-negative integer `"L"`
    - If `"L=0"`, this collapses to the White heteroskedasticity-robust
      variance estimator

If `errors = "cluster"`, then there is an additional mandatory argument:

- `cluster_var`:
  - Specifies the weighting variable
  - Must be specified if `errors = "cluster"`
  - Must be a vector
  - Must be the same length as `y`

#### Additional arguments

- `basmann`:
  - Specifies whether to use the Basmann version of the score test
  - Optional, default is `NULL`
  - Must be either `TRUE`, `FALSE` or `NULL`
- `small_sample`:
  - Specifies whether to use the small-sample variance correction
  - Optional, default is `NULL`
  - Must be either `TRUE`, `FALSE` or `NULL`
- `no_constant`:
  - Specifies whether to exclude the constant
  - Optional, default is `NULL`
  - Must be either `TRUE`, `FALSE` or `NULL`
  - Note: must be explicitly set to `"TRUE"` if constant is not wanted
    in the regression
