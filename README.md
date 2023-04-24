# Solar Fusion Reactions with Bayesian Analysis

This Python package performs analysis of nuclear reactions that are important for Solar Fusion, in particular, the associated S-factors.  The package enables the user to utilize (or add) various models for the S-factors and then will perform a Bayesian analysis of the models with the data.  The package will then perform a Bayes Model Averaging of the models utilizing the Bayes Factor as the relative weight of a model given the data.  The analysis closely follows that in Moscoso et al [arXiv:2109.00049](https://arxiv.org/abs/2109.00049), adding to that the Bayes Model Averaging.

Currently supported reactions
- $D(p,\gamma)^3{\rm He}$ -- aka $S_{12}$

See
- [Installation](#installation)
- [Usage](#usage)
- [Authors](#authors)
- [Copyright Notice](#copyright-notice)

Issues and features are tracked in the [Issues](https://github.com/walkloud/solar_fusion_reactions/issues)


### Installation
This package can be locally installed with three steps
```
git clone https://github.com/walkloud/solar_fusion_reactions
cd solar_fusion_reactions
pip install -e .
```
which should install the executable `fit_sf_reactions`.

#### Reinstalling
Reinstalling the package after updates is as simple as (from the solar_fusion_reactions directory)
```
git pull
pip install -e .
```

#### Uninstalling
To uninstall the package
```
pip uninstall solar_fusion_reactions
```

#### Testing
Test it out in a clean environment.  If you have the [Anaconda](https://www.anaconda.com/products/distribution) Python package manager installed, you can create a new env and install
```
conda create -n bare3.8 python=3.8
conda activate bare3.8
pip install -e .
```
- NOTE 1: On my new M2 mac, with a new homebrew, homebrew installs to `/opt/homebrew` rather than `/usr/local`.  Therefore, in order to get `pip` to find the HDF5 installation (after `brew install hdf5`), I need to add `export HDF5_DIR=/opt/homebrew/opt/hdf5` to my `.bash_profile`, otherwise, `pip install tables` fails.

- NOTE 2: The installation will claim to have failed, but this is a known issue with required packages described in the `setup.cfg` file and `setuptools`, not being able to find packages that where just installed as a dependency.  A test of whether it worked or not is to type

```which fit_sf_reactions```

If that returns a location, the installation should have worked.  This can be more thoroughly tested by running the help command

```fit_sf_reactions -h```

which returns the [usage](#usage) message.



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


### Copyright Notice

Low Energy Analysis of NuclEar Reactions (LEANER) Copyright (c) 2023, 
The Regents of the University of California, through Lawrence Berkeley 
National Laboratory (subject to receipt of any required approvals from 
the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.


