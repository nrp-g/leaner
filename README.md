# Solar Fusion Reactions with Bayesian Analysis

This Python package performs analysis of nuclear reactions that are important for Solar Fusion, in particular, the associated S-factors.  The package enables the user to utilize (or add) various models for the S-factors and then will perform a Bayesian analysis of the models with the data.  The package will then perform a Bayes Model Averaging of the models utilizing the Bayes Factor as the relative weight of a model given the data.  The analysis closely follows that in Moscoso et al [arXiv:2109.00049](https://arxiv.org/abs/2109.00049), adding to that the Bayes Model Averaging.

Currently supported reactions
- $D(p,\gamma)^3{\rm He}$ -- aka $S_{12}$

See
- [Installation](#installation)
- [Usage](#usage)
- [Authors](#authors)

Issues and features are tracked in the [Issues](https://github.com/walkloud/solar_fusion_reactions/issues)


### Installation
This package can be locally installed with three steps
```
git clone https://github.com/walkloud/solar_fusion_reactions
cd solar_fusion_reactions
pip install -e .
```
which should install the executable `fit_sf_reactions`.

Reinstalling the package after updates is as simple as (from the solar_fusion_reactions directory)
```
git pull
pip install -e .
```

### Usage
The analysis package uses an input file where the user specifies information such as the data file and priors for the parameters.  In addition, many options for turning on and off various options (such as adding unknown extrinsic uncertainties, the log-normal normalization factors, which models to use etc.) can be controlled at run time with the command line.  For example

```
$ fit_sf_reactions -h
usage: fit_sf_reactions [-h] [--SF_reaction SF_REACTION] [--fit_params FIT_PARAMS]
                        [--d_file D_FILE] [--d_sets D_SETS [D_SETS ...]]
                        [--models MODELS [MODELS ...]] [--f_norm]
                        [--extrinsic EXTRINSIC [EXTRINSIC ...]]
                        [--pheno_file PHENO_FILE] [--run_analysis] [--redo_fits]
                        [--report_fits] [--report_models]
                        [--prior_width PRIOR_WIDTH [PRIOR_WIDTH ...]]
                        [--prior_range PRIOR_RANGE [PRIOR_RANGE ...]] [--show_plot]
                        [--interact]

Perform Bayes analysis of nuclear reaction data S-factors

optional arguments:
  -h, --help            show this help message and exit
  --SF_reaction SF_REACTION
                        specify SF reaction [S_12]
  --fit_params FIT_PARAMS
                        specify user input file [input/fit_params.py]
  --d_file D_FILE       specify data file [data/solar_fusion_reacions.h5]
  --d_sets D_SETS [D_SETS ...]
                        user can specify d_set list to override input file
  --models MODELS [MODELS ...]
                        overide models in input file with this list
  --f_norm              S -> f S where f is an normalizationfactor to be determined
                        [True]
  --extrinsic EXTRINSIC [EXTRINSIC ...]
                        list of extrinsic statistical uncertainty modelsin analysis,
                        options are rel, abs and '' [['rel']]
  --pheno_file PHENO_FILE
                        what .dat file to use for pheno function input
                        [data/Marcucci_etal_PRL116_2016_1510.07877.dat]
  --run_analysis        run Bayes Model Analysis? [True]
  --redo_fits           redo fits even if saved? [False]
  --report_fits         print results from each model [False]
  --report_models       report model weights? [True]
  --prior_width PRIOR_WIDTH [PRIOR_WIDTH ...]
                        list priors to simultaneously perform a width study
  --prior_range PRIOR_RANGE [PRIOR_RANGE ...]
                        min, max, ds: to be used as np.arange(min,max+dm,dm) to be used
                        for prior_width study
  --show_plot           show plots? [True]
  --interact            open IPython instance after to interact with results? [False]
  ```




### Authors
This analysis package was written by André Walker-Loud (@walkloud) with consultation from G. Peter Lepage (@gplepage).
