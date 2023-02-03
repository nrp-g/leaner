import importlib
import sys
import os
import numpy as np
import gvar as gv
import lsqfit
import scipy.stats as stats
import matplotlib.pyplot as plt

import fitter.load_data as ld
import fitter.plotting as plot
import fitter.sf_fit_functions as sf_fit

class FitManager():

    def __init__(self,args):
        # load the data
        self.x, self.y  = ld.load_data(args.d_file, args.d_path, args.d_sets)
        self.f_norm = args.f_norm
        self.args = args

        # load the fit_params from user input file
        sys.path.append(os.path.dirname(os.path.abspath(args.fit_params)))
        fp = importlib.import_module(args.fit_params.split('/')[-1].split('.py')[0])
        self.fit_params = fp.fit_params[args.SF_reaction]

    def plot_data(self):
        self.ax = plot.plot_data(self.x, self.y, self.fit_params['plot_params'])

    def fit_extrinsic_sig(self,z):
        self.new_y = {}
        plausibility = 0
        xcutoff = 1.
        #xcutoff = 10.
        for k in z:
            #dy = gv.mean(self.y[k]) * z[k]**2
            dy = xcutoff * np.arctan(z[k]**2 / xcutoff) * gv.mean(self.y[k])
            #dy = z[k]**2
            #np.array([d.mean for d in self.y[k]]) * z[k]
            self.new_y[k] = self.y[k] + gv.gvar(np.zeros_like(dy), dy)
            plausibility -= z[k]**2 / (2*xcutoff)
        fit_dict = {
            'data' : (self.x, self.new_y),
            'fcn'  : self.tmp_model_fit.fit_func,
            'p0'   : self.tmp_p0,
            'prior': self.tmp_p
        }
        return fit_dict
        #return (fit_dict, plausibility)

    def fit_models(self):
        self.fit_results = {}
        for model in self.fit_params['models']:
            self.tmp_model_fit = sf_fit.SFFunctions(model, f_norm=self.f_norm)
            self.tmp_p,self.tmp_p0 = self.tmp_model_fit.prune_priors(self.fit_params['priors'])
            if self.args.extrinsic:
                print(model, 'extrinsic hunt')
                ffit,z = lsqfit.empbayes_fit({'Turkat_2021':0.1, 'Mossa_2020':0.1,
                            'Tisma_2019':0.1, 'Casella_2002':0.1, 'Schmid_1997':0.1,
                            'Ma_1997':0.1,'Warren_1963':0.1}, self.fit_extrinsic_sig)
                print(z)
                self.fit_results[model] = ffit
            else:
                self.fit_results[model] = lsqfit.nonlinear_fit(
                                            data  = (self.x, self.y),
                                            prior = self.tmp_p,
                                            p0    = self.tmp_p0,
                                            fcn   = self.tmp_model_fit.fit_func,
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
        self.model_var = 0
        for model in self.fit_params['models']:
            model_fit = sf_fit.SFFunctions(model, f_norm=False)
            tmp = model_fit.fit_func(self.x_plot, self.fit_results[model].p)
            mean = np.array([k.mean for k in tmp['plot']])
            sdev = np.array([k.sdev for k in tmp['plot']])

            self.model_avg += gv.gvar(self.weight[model] * mean,
                                      np.sqrt(self.weight[model]) * sdev)
            self.model_var += self.weight[model] * mean**2
        self.model_var += -np.array([k.mean**2 for k in self.model_avg])


    def plot_fit(self):
        self.do_model_avg()
        yy   = np.array([k.mean for k in self.model_avg])
        var  = np.array([k.var  for k in self.model_avg])
        var += self.model_var
        dy   = np.sqrt(var)
        self.ax.fill_between(self.x_plot['plot'], yy-dy, yy+dy, color='k', alpha=.3)

    def report_fits(self):
        for model in self.fit_params['models']:
            print('-----------------------------------------------------------')
            print(model)
            print('-----------------------------------------------------------')
            print(self.fit_results[model])

    def plot_pdf(self):
        pdf = 0.
        cdf = 0.
        pdf_x = np.arange(1.5e-7,2.4e-7,1e-9)
        self.get_weights()
        for model in self.fit_params['models']:
            w_i = self.weight[model]
            S_0 = self.fit_results[model].p['S_0']
            p   = stats.norm.pdf(pdf_x, S_0.mean, S_0.sdev)
            pdf += w_i * p
            cdf += w_i * stats.norm.cdf(pdf_x, S_0.mean, S_0.sdev)
        self.pdf = pdf
        self.cdf = cdf
        hist = plt.figure('hist')
        ax   = plt.axes([0.12, 0.12, 0.85, 0.85])
        ax.plot(pdf_x, self.pdf, color='k')
        ax.fill_between(x=pdf_x, y1=self.pdf, color='k', alpha=.2)
        ax.set_xlabel(r'$S_0$[MeV b]', fontsize=16)
        ax.set_ylabel(r'likelihood', fontsize=16)
        if not os.path.exists('figures'):
            os.makedirs('figures')
        plt.savefig('figures/S_0_hist.pdf',transparent=True)
