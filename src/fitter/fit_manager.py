import importlib
import sys
import os
import numpy as np
import gvar as gv
import lsqfit

import fitter.load_data as ld
import fitter.plotting as plot
import fitter.sf_fit_functions as sf_fit

class FitManager():

    def __init__(self,args):
        # load the data
        self.x, self.y  = ld.load_data(args.d_file, args.d_path, args.d_sets)

        # load the fit_params from user input file
        sys.path.append(os.path.dirname(os.path.abspath(args.fit_params)))
        fp = importlib.import_module(args.fit_params.split('/')[-1].split('.py')[0])
        self.fit_params = fp.fit_params[args.SF_reaction]

    def plot_data(self):
        self.ax = plot.plot_data(self.x, self.y, self.fit_params['plot_params'])

    def fit_models(self):
        self.fit_results = {}
        for model in self.fit_params['models']:
            model_fit = sf_fit.SFFunctions(model)
            p,p0 = model_fit.prune_priors(self.fit_params['priors'])
            self.fit_results[model] = lsqfit.nonlinear_fit(
                                    data  = (self.x, self.y),
                                    prior = p,
                                    p0    = p0,
                                    fcn   = model_fit.fit_func,
                                    )

    def plot_fit(self):
        for model in self.fit_params['models']:
            model_fit = sf_fit.SFFunctions(model)
            x = {}
            xi = self.fit_params['plot_params']['x_lim'][0]
            xf = self.fit_params['plot_params']['x_lim'][1]
            x['plot'] = np.arange(xi, xf, xi/10)
            result = model_fit.fit_func(x, self.fit_results[model].p)
            yy = np.array([k.mean for k in result['plot']])
            dy = np.array([k.sdev for k in result['plot']])
            self.ax.fill_between(x['plot'], yy-dy, yy+dy, color='k', alpha=.3)
