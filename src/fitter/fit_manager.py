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

    def get_weights(self):
        logGBF = {}
        for model in self.fit_params['models']:
            logGBF[model] = self.fit_results[model].logGBF
        logGBF_max  = max([logGBF[k] for k in logGBF])
        d_logGBF    = {k:logGBF[k]-logGBF_max for k in logGBF}
        norm        = np.sum([np.exp(d_logGBF[k]) for k in logGBF])
        self.weight = {k:np.exp(d_logGBF[k])/norm for k in logGBF}


    def do_model_avg(self):
        self.x_plot = {}
        xi = self.fit_params['plot_params']['x_lim'][0]
        xf = self.fit_params['plot_params']['x_lim'][1]
        self.x_plot['plot'] = np.arange(xi, xf, xi/10)
        self.get_weights()
        self.model_avg = 0
        for model in self.fit_params['models']:
            model_fit = sf_fit.SFFunctions(model, norm=False)
            tmp = model_fit.fit_func(self.x_plot, self.fit_results[model].p)
            self.model_avg += self.weight[model] * tmp['plot']


    def plot_fit(self):
        self.do_model_avg()
        yy = np.array([k.mean for k in self.model_avg])
        dy = np.array([k.sdev for k in self.model_avg])
        self.ax.fill_between(self.x_plot['plot'], yy-dy, yy+dy, color='k', alpha=.3)

    def report_fits(self):
        for model in self.fit_params['models']:
            print('-----------------------------------------------------------')
            print(model)
            print('-----------------------------------------------------------')
            print(self.fit_results[model])
