# oidrobust
The `oidrobust` package provides facilities for robust overidentification testing in linear IV models, currently available in R and python. Also, MATLAB files are available to replicate the simulations and empirical applications:

Lane, S., & Windmeijer, F. (2025). Overidentification testing with weak instruments and heteroskedasticity. arXiv preprint arXiv:2509.21096.

## Installation

### R

Install the R development version from GitHub:
```r
# install.packages("devtools")
devtools::install_github("stuart-lane/oidrobust", subdir="R")
```
### Python
Install the python development version from GitHub:
```python
pip install git+https://github.com/stuart-lane/oidrobust.git#subdirectory=Python
```

 ## ⚠️ Development Status
 
This package is in the late stages of development, but has not undergone extensive testing. While I believe the implementation is correct, bugs may exist. Use with caution and please report any issues.
