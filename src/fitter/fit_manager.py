import importlib
import sys
import os
import numpy as np
import gvar as gv
import lsqfit
import scipy.stats as stats
import matplotlib.pyplot as plt
from glob import glob

import fitter.load_data as ld
import fitter.plotting as plot
import fitter.sf_fit_functions as sf_fit

class FitManager():

    def __init__(self,args):
        # load the fit_params from user input file
        sys.path.append(os.path.dirname(os.path.abspath(args.fit_params)))
        fp = importlib.import_module(args.fit_params.split('/')[-1].split('.py')[0])
        self.fit_params = fp.fit_params[args.SF_reaction]

        # load the data
        if args.d_sets:
            dsets = args.d_sets
        else:
            dsets = self.fit_params['dsets']
        self.x, self.y  = ld.load_data(args.d_file, args.d_path, dsets)
        self.f_norm = args.f_norm
        self.args = args


    def plot_data(self):
        self.ax = plot.plot_data(self.x, self.y, self.fit_params['plot_params'])

    def fit_extrinsic_sig(self,z):
        self.new_y = {}
        plausibility = 0
        xcutoff = 1.
        #xcutoff = 10.
        for k in z:
            #dy = gv.mean(self.y[k]) * z[k]**2
            #dy = xcutoff * np.arctan(z[k]**2 / xcutoff)
            dy = xcutoff * np.arctan(z[k]**2 / xcutoff) * gv.mean(self.y[k])
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
            self.tmp_p, self.tmp_p0 = self.tmp_model_fit.prune_priors(self.fit_params['priors'],self.y)

            if self.args.extrinsic:
                print(model, 'extrinsic hunt')
                z0 = {k:v for k,v in self.fit_params['z0'].items() if k in self.y}
                ffit,z = lsqfit.empbayes_fit(z0, self.fit_extrinsic_sig)
                print(z)
                self.fit_results[model] = ffit
            else:
                self.fit_results[model] = lsqfit.nonlinear_fit(
                                            data  = (self.x, self.y),
                                            prior = self.tmp_p,
                                            p0    = self.tmp_p0,
                                            fcn   = self.tmp_model_fit.fit_func,
                                            )

    def prior_width(self):
        self.prior_width_results = {}
        for model in self.fit_params['models']:
            self.tmp_model_fit = sf_fit.SFFunctions(model, f_norm=self.f_norm)
            self.tmp_p, self.tmp_p0 = self.tmp_model_fit.prune_priors(self.fit_params['priors'],self.y)
            pi, pf, dp = self.args.prior_range
            for sp in np.arange(pi, pf+dp, dp):
                model_save = model.replace(' ','_')
                for p in self.tmp_p:
                    if p in self.args.prior_width:
                        self.tmp_p[p] = gv.gvar(self.tmp_p[p].mean, sp)
                        model_save += '_'+p
                model_save += '_dS_%.2e' %(sp)
                print(model_save)
                if not os.path.exists('pickled_fits'):
                    os.makedirs('pickled_fits')
                if not os.path.exists('pickled_fits/'+model_save+'.p'):
                    if self.args.extrinsic:
                        z0 = {k:v for k,v in self.fit_params['z0'].items() if k in self.y}
                        ffit,z = lsqfit.empbayes_fit(z0, self.fit_extrinsic_sig)
                        self.prior_width_results[model_save] = ffit
                    else:
                        self.prior_width_results[model_save] = lsqfit.nonlinear_fit(
                                                    data  = (self.x, self.y),
                                                    prior = self.tmp_p,
                                                    p0    = self.tmp_p0,
                                                    fcn   = self.tmp_model_fit.fit_func,
                                                    )
                    gv.dump(self.prior_width_results[model_save], 'pickled_fits/'+model_save+'.p', add_dependencies=True)
                else:
                    print(model_save+' already exists - load fit')
                    self.prior_width_results[model_save] = gv.load('pickled_fits/'+model_save+'.p')

            fit_search = model_save.split('_dS')[0]
            all_fits   = glob('pickled_fits/'+fit_search+'*.p')
            for fit_result in all_fits:
                fit = fit_result.split('/')[-1].split('.p')[0]
                if fit not in self.prior_width_results:
                    self.prior_width_results[fit] = gv.load(fit_result)
            logGBF = {}
            ds = {}
            for fit in self.prior_width_results:
                logGBF[fit] = self.prior_width_results[fit].logGBF
                ds[fit]     = float(fit.split('dS_')[-1])
            logGBF_max = max([logGBF[k] for k in logGBF])
            d_logGBF   = {k:logGBF[k]-logGBF_max for k in logGBF}
            norm       = np.sum([np.exp(d_logGBF[k]) for k in logGBF])
            weight     = {k:np.exp(d_logGBF[k])/norm for k in logGBF}
            fits       = [fit for fit in weight]
            x          = np.array([ds[fit] for fit in fits])
            i_x        = np.argsort(x)
            w          = np.array([weight[fit] for fit in fits])
            plt.figure(model+' prior width')
            ax = plt.axes([.12, .12, .87, .87])
            ax.plot(x[i_x],w[i_x],marker='s')
            ax.text(.05, .9, model, transform=ax.transAxes,fontsize=20)
            ax.set_xlabel(r'$\tilde{\sigma}$', fontsize=16)
            ax.set_ylabel('relative weight', fontsize=16)
            plt.savefig('figures/'+model+'_prior_width_study.pdf', transparent=True)

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


    def model_avg_S(self,E,plot_hist=False):
        self.get_weights()
        E_result = {}
        E_result['eval'] = np.array([E])
        mean = 0
        var  = 0
        for model in self.fit_params['models']:
            model_fit = sf_fit.SFFunctions(model, f_norm=False)
            tmp   = model_fit.fit_func(E_result, self.fit_results[model].p)
            t_m   = np.array([k.mean for k in tmp['eval']])
            t_s   = np.array([k.sdev for k in tmp['eval']])
            mean += gv.gvar(self.weight[model] * t_m, np.sqrt(self.weight[model]) * t_s)
            var  += self.weight[model] * t_m**2
        var += -np.array([k.mean**2 for k in mean])
        total_var = mean[0].var + var
        print('S(%f) = %s +- %.1e' %(E,mean[0],np.sqrt(var)[0]))

        if plot_hist:
            pdf = 0.
            cdf = 0.
            pdf_x = np.arange(mean[0].mean-10*np.sqrt(total_var), 
                              mean[0].mean+10*np.sqrt(total_var),
                              np.sqrt(total_var)/33)
            self.pdf_model = {}
            for model in self.fit_params['models']:
                w_i       = self.weight[model]
                model_fit = sf_fit.SFFunctions(model, f_norm=False)
                S_model   = model_fit.fit_func(E_result, self.fit_results[model].p)
                p         = stats.norm.pdf(pdf_x, S_model['eval'][0].mean, S_model['eval'][0].sdev)
                pdf      += w_i * p
                self.pdf_model[model] = w_i * p
                cdf      += w_i * stats.norm.cdf(pdf_x, S_model['eval'][0].mean, S_model['eval'][0].sdev)
            if 'pdf' not in dir(self):
                self.pdf = {}
                self.cdf = {}
            self.pdf[E] = pdf
            self.cdf[E] = cdf
            hist = plt.figure('hist E=%f' %E)
            ax   = plt.axes([0.12, 0.12, 0.85, 0.85])
            ax.plot(pdf_x, self.pdf[E], color='k')
            ax.fill_between(x=pdf_x, y1=self.pdf[E], color='k', alpha=.2)
            for model in self.fit_params['models']:
                ax.fill_between(pdf_x, self.pdf_model[model], alpha=.4, label=model)
            ax.set_xlabel(r'$dS(E=%.4f)$[MeV b]' %(E), fontsize=16)
            ax.set_ylabel(r'$d\mathcal{P}$', fontsize=16)
            ax.legend()
            if not os.path.exists('figures'):
                os.makedirs('figures')
            plt.savefig('figures/S_E%.4f_hist.pdf' %(E),transparent=True)

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
        pdf_x = np.arange(1.1e-7,3.0e-7,1e-9)
        self.get_weights()
        self.pdf_model = {}
        for model in self.fit_params['models']:
            w_i = self.weight[model]
            S_0 = self.fit_results[model].p['S_0']
            p   = stats.norm.pdf(pdf_x, S_0.mean, S_0.sdev)
            pdf += w_i * p
            self.pdf_model[model] = w_i * p
            cdf += w_i * stats.norm.cdf(pdf_x, S_0.mean, S_0.sdev)
        self.pdf = pdf
        self.cdf = cdf
        hist = plt.figure('hist')
        ax   = plt.axes([0.12, 0.12, 0.85, 0.85])
        ax.plot(pdf_x, self.pdf, color='k')
        ax.fill_between(x=pdf_x, y1=self.pdf, color='k', alpha=.2)
        for model in self.fit_params['models']:
            ax.fill_between(pdf_x, self.pdf_model[model], alpha=.4, label=model)
        ax.set_xlabel(r'$dS_0$[MeV b]', fontsize=16)
        ax.set_ylabel(r'$d\mathcal{P}$', fontsize=16)
        ax.legend()
        if not os.path.exists('figures'):
            os.makedirs('figures')
        plt.savefig('figures/S_0_hist.pdf',transparent=True)
