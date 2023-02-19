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
The analysis package uses an input file where the user specifies information such as the data file and priors for the parameters.




### Authors
This analysis package was written by André Walker-Loud (@walkloud) with consultation from G. Peter Lepage (@gplepage).
