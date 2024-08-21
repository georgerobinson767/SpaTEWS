# SpaTEWS package

This MATLAB package computes early warning signals (EWS) for spatial/multivariate and temporal data sets, with the primary purpose of comparing temporal and spatially informed EWS for spatial or multivariate data sets. This package has been designed for use on both empirical and synthetic data sets. We provide the synthetic_data.m file to generate synthetic data sets to explore using the functions available in this package.

## Functions
* **run_ews.m** - computes EWS and/or performs sensitivity analysis
* sens_fig.m - produces sensitivity figures
* ews_fig.m - produces EWS figures
* ews_calculator.m - computes EWS
* synthetic_data.m - produces synthetic data sets

Also included in this package is a slightly edited version of
* Modified_MannKendall_test.m

written by Atharva Aalok,
https://github.com/atharvaaalok/Modified-MannKendall-Test, which implements the Modified Mann-Kendall test proposed by Hamed and Rao (1998). Furthermore, the functions:
* run_all.m
* create_all_figs.m

make use of the functions described above and can be used to reproduce the data and some figures displayed in our paper.

## User guide
See the **demo.m** file for some examples of different ways to use the available functions in this package. Information regarding the different inputs for each function is available in both the **run_ews.m** and **synthetic_data.m** files. The documentation for using the **run_all.m** and **create_all_figs.m** files can be found at the end of the demo.m file.

## Dependencies

* Econometrics toolbox - required for Modified_MannKendall_test.m

The Econometrics toolbox by MathWorks requires
* Optimization Toolbox
* Statistics and Machine Learning Toolbox

### Optional
* Parallel Computing Toolbox (recommended)
* tight_subplot(Nh, Nw, gap, marg_h, marg_w)
* export_fig