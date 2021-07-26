# ASRE-Compute

From a pandas dataframe, this program computes the Absolute Rule Effect or cognitive biais Absolute Stochastic Rule Effect through A-learning and provides asymptotic 95% confidence intervals from M-estimation's sandwich formula.

## Instructions

1. Install:

```
pip install asre-compute
```

2. Import, specify the models, plot:

```python
from asre_compute import tools

tools.asre_package(observational_dataframe, 
          rule = "rule",
          ttt = "cabg",
          y = "Y",
          ps_predictors = ["age", "crcl_log", "copd", "tvd", "lmcad", "both"],
          pronostic_predictors = ["tvd", "lmcad", "both", "syntax", "age", "crcl", "diabetes", "insulin", "lvef", "smoking", "pvd", "copd"],
          ctst_vrb = ['syntax', 'tvd', 'lmcad'],
          est='ASRE_cb', alpha = .5, n_alphas=5, precision=3)
```

3. Bob's your uncle
