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
        self.args       = args
        if args.models:
            self.fit_params['models'] = args.models
        self.fit_params['plot_params']['residuals'] = args.residuals
        if args.save_fits:
            if not os.path.exists('pickled_fits'):
                os.makedirs('pickled_fits')


    def plot_data(self):
        if self.fit_params['plot_params']['residuals']:
            self.ax, self.ax_resid, self.ax_shift, self.ax_resid_shift =\
                plot.plot_data(self.x, self.y, self.fit_params['plot_params'],
                               error=True, marker='o', label=True, alpha=0.5)
        else:
            self.ax = plot.plot_data(self.x, self.y, self.fit_params['plot_params'],
                                     error=False, marker='s', label=False)

    def fit_extrinsic_sig_rel(self,z):
        ''' add an unknown statistical uncertainty to the data relative to the data.mean

            new_y = y + (0, y.mean * z_cutoff * arctan(z**2 / z_cutoff) )

            The parameter z will be optimized with lsqfit.empbayes_fit
        '''
        self.new_y = {}
        z_cutoff = 1.
        for k in z:
            dy = z_cutoff * np.arctan(z[k]**2 / z_cutoff) * gv.mean(self.y[k])
            self.new_y[k] = self.y[k] + gv.gvar(np.zeros_like(dy), dy)
        fit_dict = {
            'data' : (self.x, self.new_y),
            'fcn'  : self.tmp_model_fit.fit_func,
            'p0'   : self.tmp_p0,
            'prior': self.tmp_p
        }
        return fit_dict

    def fit_extrinsic_sig_abs(self,z):
        ''' add an unknown statistical uncertainty to the data absolute to the data.mean

            new_y = y + z_cutoff * arctan(z**2 / z_cutoff)

            The parameter z will be optimized with lsqfit.empbayes_fit
        '''
        self.new_y = {}
        z_cutoff = 1.
        for k in z:
            dy = z_cutoff * np.arctan(z[k]**2 / z_cutoff)
            self.new_y[k] = self.y[k] + gv.gvar(np.zeros_like(dy), dy)
        fit_dict = {
            'data' : (self.x, self.new_y),
            'fcn'  : self.tmp_model_fit.fit_func,
            'p0'   : self.tmp_p0,
            'prior': self.tmp_p
        }
        return fit_dict

    def extrinsic_sig(self, z, ext='rel'):
        ''' add extra "unknown" statistical uncertainty to data
            either as "relative"
                dy = y.mean * z_cutoff * arctan(z**2 / z_cutoff)
            or "absolute"
                dy = z_cutoff * arctan(z**2 / z_cutoff)
        '''
        dy_ext = {}
        z_cutoff = 1.
        for k in z:
            if ext == 'rel':
                dy = z_cutoff * np.arctan(z[k]**2 / z_cutoff) * gv.mean(self.y[k])
            elif ext == 'abs':
                dy = z_cutoff * np.arctan(z[k]**2 / z_cutoff) * np.ones_like(self.y[k])
            else:
                sys.exit('I do not understand your choice of extrinsic uncertainty,', ext)
            dy_ext[k] = gv.gvar(np.zeros_like(dy), dy)
        
        return dy_ext

    def fit_models(self):
        self.fit_results = {}
        for model in self.fit_params['models']:
            for extrinsic in self.args.extrinsic:
                if extrinsic:
                    full_m = model+'_extrinsic_sig_'+extrinsic
                else:
                    full_m = model+'_extrinsic_sig_none'

                offset = 'offset' in model
                self.tmp_model_fit = sf_fit.SFFunctions(model, f_norm=self.args.f_norm, 
                                                        offset=offset, pheno_file=self.args.pheno_file)
                self.tmp_p, self.tmp_p0 = self.tmp_model_fit.prune_priors(self.fit_params['priors'], self.y)

                # if the fit is saved, load it, otherwise, do the analysis
                saved_fit = 'pickled_fits/'+self.args.SF_reaction+'_'+full_m+'.p'

                if self.args.redo_fits and os.path.exists(saved_fit):
                    os.remove(saved_fit)
                if os.path.exists(saved_fit):
                    print('fit already performed - loading')
                    self.fit_results[full_m] = gv.load(saved_fit)
                else:
                    if extrinsic:
                        print(model, 'optimize extrinsic uncertainty')
                        z0 = {k:v for k,v in self.fit_params['z0'].items() if k in self.y}
                        ext_model = 'fit_extrinsic_sig_'+extrinsic

                        ffit,z = lsqfit.empbayes_fit(z0, getattr(self, ext_model))
                        
                        ffit.extrinsic_sig = z
                        self.fit_results[full_m] = ffit
                    else:
                        self.fit_results[full_m] = lsqfit.nonlinear_fit(
                                                    data  = (self.x, self.y),
                                                    prior = self.tmp_p,
                                                    p0    = self.tmp_p0,
                                                    fcn   = self.tmp_model_fit.fit_func,
                                                    )
                    if self.args.save_fits:
                        gv.dump(self.fit_results[full_m], saved_fit, add_dependencies=True)

    def prior_width(self):
        self.prior_width_results = {}
        for model in self.fit_params['models']:
            for extrinsic in self.args.extrinsic:
                if extrinsic:
                    full_m = model+'_extrinsic_sig_'+extrinsic
                else:
                    full_m = model+'_extrinsic_sig_none'

                offset = 'offset' in model
                self.tmp_model_fit = sf_fit.SFFunctions(model, f_norm=self.args.f_norm, 
                                                        offset=offset, pheno_file=self.args.pheno_file)
                self.tmp_p, self.tmp_p0 = self.tmp_model_fit.prune_priors(self.fit_params['priors'],self.y)

                pi, pf, dp = self.args.prior_range
                for sp in np.arange(pi, pf+dp, dp):
                    model_key = str(full_m)
                    for p in self.tmp_p:
                        if p in self.args.prior_width:
                            self.tmp_p[p] = gv.gvar(self.tmp_p[p].mean, sp)
                            model_key += '_'+p
                    model_key += '_dS_%.2e' %(sp)
                    model_save = 'pickled_fits/'+self.args.SF_reaction+'_'+model_key+'.p'
                    print(model_key)
                    if not os.path.exists('pickled_fits'):
                        os.makedirs('pickled_fits')
                    if not os.path.exists(model_save):
                        if extrinsic:
                            z0 = {k:v for k,v in self.fit_params['z0'].items() if k in self.y}
                            ext_model = 'fit_extrinsic_sig_'+extrinsic
                            ffit,z = lsqfit.empbayes_fit(z0, getattr(self, ext_model))
                            ffit.extrinsic_siz = z
                            self.prior_width_results[model_key] = ffit
                        else:
                            self.prior_width_results[model_key] = lsqfit.nonlinear_fit(
                                                        data  = (self.x, self.y),
                                                        prior = self.tmp_p,
                                                        p0    = self.tmp_p0,
                                                        fcn   = self.tmp_model_fit.fit_func,
                                                        )
                        if self.args.save_fits:
                            gv.dump(self.prior_width_results[model_key], model_save, add_dependencies=True)
                    else:
                        print(model_key+' already exists - load fit')
                        self.prior_width_results[model_key] = gv.load(model_save)

                all_fits = glob('pickled_fits/'+full_m+'*.p')
                for fit_result in all_fits:
                    fit = fit_result.split('/')[-1].split('.p')[0]
                    if full_m in fit and fit not in self.prior_width_results:
                        self.prior_width_results[fit] = gv.load(fit_result)
                logGBF = {}
                ds = {}
                for fit in self.prior_width_results:
                    if full_m in fit:
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
                plt.figure(full_m+' prior width')
                ax = plt.axes([.12, .12, .87, .87])
                ax.plot(x[i_x],w[i_x],marker='s')
                ax.set_ylim(-0.005, 1.15*w.max())
                ax.text(.05, .9, full_m, transform=ax.transAxes,fontsize=20)
                ax.set_xlabel(r'$\tilde{\sigma}$', fontsize=16)
                ax.set_ylabel('relative weight', fontsize=16)
                plt.savefig('figures/'+full_m+'_prior_width_study.pdf', transparent=True)

    def get_weights(self):
        logGBF = {}
        for result in self.fit_results:
            model = result.split('_')[0]+'_'+result.split('_')[1]
            logGBF[result] = self.fit_results[result].logGBF
        logGBF_max  = max([logGBF[k] for k in logGBF])
        d_logGBF    = {k:logGBF[k]-logGBF_max for k in logGBF}
        norm        = np.sum([np.exp(d_logGBF[k]) for k in logGBF])
        self.weight = {k:np.exp(d_logGBF[k])/norm for k in logGBF}


    def do_model_avg(self):
        self.x_plot = {}
        xi = self.fit_params['plot_params']['x_lim'][0]
        if not any(['pheno' in model for model in self.fit_results]):
            xf = self.fit_params['plot_params']['x_lim'][1]
        else:
            xf = 2
        self.x_plot['plot'] = np.arange(xi, xf, xi/10)
        self.get_weights()
        self.model_avg = 0
        self.model_var = 0
        for result in self.fit_results:
            model     = result.split('_')[0]+'_'+result.split('_')[1]
            offset    = 'offset' in model
            model_fit = sf_fit.SFFunctions(model, f_norm=False, 
                                           offset=offset, pheno_file=self.args.pheno_file)
            tmp       = model_fit.fit_func(self.x_plot, self.fit_results[result].p)
            mean      = np.array([k.mean for k in tmp['plot']])
            sdev      = np.array([k.sdev for k in tmp['plot']])

            self.model_avg += gv.gvar(self.weight[result] * mean,
                                      np.sqrt(self.weight[result]) * sdev)
            self.model_var += self.weight[result] * mean**2
        self.model_var += -np.array([k.mean**2 for k in self.model_avg])

    def report_models(self):
        self.get_weights()
        models  = [model for model in self.fit_results]
        weights = [self.weight[model] for model in models]
        i_w     = np.argsort(weights)[::-1]
        print('')
        print('BAYES MODEL AVERAGE: rel. weight = exp(logGBF_i - logGBF_max)')
        print('----------------------------------------------------------------------------------------------')
        print('%30s  %6s  %8s  chi^2/dof [dof]  Q' %('model', 'logGBF', ' weight '))
        print('----------------------------------------------------------------------------------------------')
        for m in i_w:
            model  = models[m]
            logGBF = self.fit_results[model].logGBF
            chi2   = self.fit_results[model].chi2
            dof    = self.fit_results[model].dof
            Q      = self.fit_results[model].Q
            w      = weights[m]
            if Q >= 0.001:
                print("%30s  %.1f  %.2e    %.4f  [%d]   %.3f" %(model, logGBF, w, chi2/dof, dof, Q))
            else:
                print("%30s  %.1f  %.2e    %.4f  [%d]   %.2e" %(model, logGBF, w, chi2/dof, dof, Q))
        print('----------------------------------------------------------------------------------------------\n')

    def model_avg_S(self, E, print_result=True, plot_hist=False):
        self.get_weights()
        E_result = {}
        E_result['eval'] = np.array([E])
        mean = 0
        var  = 0
        for result in self.fit_results:
            model     = result.split('_')[0]+'_'+result.split('_')[1]
            offset    = 'offset' in model
            model_fit = sf_fit.SFFunctions(model, f_norm=False, 
                                           offset=offset, pheno_file=self.args.pheno_file)
            tmp       = model_fit.fit_func(E_result, self.fit_results[result].p)
            t_m       = np.array([k.mean for k in tmp['eval']])
            t_s       = np.array([k.sdev for k in tmp['eval']])
            mean     += gv.gvar(self.weight[result] * t_m, np.sqrt(self.weight[result]) * t_s)
            var      += self.weight[result] * t_m**2
        var += -np.array([k.mean**2 for k in mean])
        total_var = mean[0].var + var
        self.model_var = var
        self.mean_var  = mean[0].var

        S_ma = gv.gvar(mean[0].mean, np.sqrt(total_var)[0])
        if print_result:
            print('S(E=%f) = %s +- %.1e' %(E,mean[0],np.sqrt(var)[0]))

        if plot_hist:
            pdf = 0.
            cdf = 0.
            pdf_x = np.arange(mean[0].mean-10*np.sqrt(total_var), 
                              mean[0].mean+10*np.sqrt(total_var),
                              np.sqrt(total_var)/33)
            self.pdf_model = {}
            for result in self.fit_results:
                model     = result.split('_')[0]+'_'+result.split('_')[1]
                offset    = 'offset' in model
                w_i       = self.weight[result]
                model_fit = sf_fit.SFFunctions(model, f_norm=False, 
                                               offset=offset, pheno_file=self.args.pheno_file)
                S_model   = model_fit.fit_func(E_result, self.fit_results[result].p)
                p         = stats.norm.pdf(pdf_x, S_model['eval'][0].mean, S_model['eval'][0].sdev)
                pdf      += w_i * p
                self.pdf_model[result] = w_i * p
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
            for result in self.fit_results:
                ax.fill_between(pdf_x, self.pdf_model[result], alpha=.4, label=result)
            ax.set_xlabel(r'$dS(E=%.4f)$[MeV b]' %(E), fontsize=16)
            ax.set_ylabel(r'$d\mathcal{P}$', fontsize=16)
            ax.legend()
            if not os.path.exists('figures'):
                os.makedirs('figures')
            plt.savefig('figures/S_E%.4f_hist.pdf' %(E),transparent=True)

        return S_ma

    def ma_extrinsic(self):
        self.dy_extrinsic_ma = gv.BufferDict()
        for dset in self.x:
            self.dy_extrinsic_ma[dset] = np.zeros_like(self.y[dset])
            for model in self.fit_results:
                ext = model.split('_')[-1]
                dy = self.extrinsic_sig(self.fit_results[model].extrinsic_sig, ext=ext)
                self.dy_extrinsic_ma[dset] += self.weight[model] * dy[dset]
        self.y_extrinsic_ma = gv.BufferDict()
        for k in self.y:
            self.y_extrinsic_ma[k] = self.y[k] + self.dy_extrinsic_ma[k]

    def ma_f_norm(self):
        self.f_inv = gv.BufferDict()
        for dset in self.x:
            for model in self.fit_results:
                if dset in self.f_inv:
                    self.f_inv[dset] += self.weight[model] / self.fit_results[model].p['f_'+dset]
                else:
                    self.f_inv[dset]  = self.weight[model] / self.fit_results[model].p['f_'+dset]

    def plot_shifted_y(self):
        for dset in self.x:
            x = self.x[dset]
            clr = self.fit_params['plot_params']['colors'][dset]
            lbl = dset.replace('_','\_')

            extrinsic=True
            if extrinsic and self.args.f_norm:
                # full uncertainty
                y_new = self.y_extrinsic_ma[dset] * self.f_inv[dset]
                y_m = [k.mean for k in y_new]
                y_s = [k.sdev for k in y_new]
                self.ax_shift.errorbar(x, y_m, yerr=y_s,
                                 marker='s', mfc='None', c=clr, alpha=0.5, linestyle='None')
                
                '''
                # only extrinsic uncertainty
                y_new = self.y_extrinsic_ma[dset] * self.f_inv[dset].mean
                y_m = [k.mean for k in y_new]
                y_s = [k.sdev for k in y_new]
                self.ax.errorbar(x, y_m, yerr=y_s,
                                 marker='None', c='k', alpha=0.5, linestyle='None')
                
                # only original stoch uncertainty
                y_new = self.y[dset] * self.f_inv[dset].mean
                y_m = [k.mean for k in y_new]
                y_s = [k.sdev for k in y_new]
                self.ax.errorbar(x, y_m, yerr=y_s,
                                 marker='o', c=clr, mfc='None', label=lbl, linestyle='None')
                '''
        if 'legend_loc' in self.fit_params['plot_params']:
            self.ax.legend(loc=self.fit_params['plot_params']['legend_loc'])
        else:
            self.ax.legend()

    def plot_fit(self):
        self.do_model_avg()
        yy   = np.array([k.mean for k in self.model_avg])
        var  = np.array([k.var  for k in self.model_avg])
        var += self.model_var
        dy   = np.sqrt(var)
        self.ax.fill_between(self.x_plot['plot'], yy-dy, yy+dy, color='k', alpha=.3)
        self.ax_shift.fill_between(self.x_plot['plot'], yy-dy, yy+dy, color='k', alpha=.3)

        if self.fit_params['plot_params']['residuals']:
            # mean.var + model_var
            var = self.model_var + self.mean_var
            dy  = np.sqrt(var) #''' divide by mean '''
            self.ax_resid.fill_between(self.x_plot['plot'], -dy, dy, color='k', alpha=.2)
            # mean.var
            var = self.mean_var
            dy  = np.sqrt(var)
            self.ax_resid.fill_between(self.x_plot['plot'], -dy, dy, color='k', alpha=.4)
            self.ax_resid_shift.fill_between(self.x_plot['plot'], -dy, dy, color='k', alpha=.4)

            # subtract data from fit
            plot_params = self.fit_params['plot_params']

            # plot shifted data
            # determine data with extrinsic uncertainty
            if any([k in ['abs','rel'] for k in self.args.extrinsic]):
                extrinsic = True
                self.ma_extrinsic()
            else:
                extrinsic = False
            # determine model avg normalization
            if self.args.f_norm:
                self.ma_f_norm()

            for dset in self.x:
                x = self.x[dset]
                clr  = plot_params['colors'][dset]
                for i_e, E in enumerate(x):
                    y_ma = self.model_avg_S(E, print_result=False)

                    dy = (self.y[dset][i_e] - y_ma ) / y_ma
                    if not (extrinsic or self.args.f_norm):
                        self.ax_resid.errorbar(E, dy.mean, yerr=dy.sdev,
                                               marker='o', c=clr, mfc='None')
                    else:
                        # plot original data on left
                        self.ax_resid.errorbar(E, dy.mean, yerr=dy.sdev,
                                               marker='o', c=clr, mfc='None')

                        # plot original data without error bar
                        #self.ax_resid.plot(E, dy.mean, alpha=0.4,
                        #                   marker='s', c=clr, mfc='None')

                        # if we have extrinsic and f_norm
                        if extrinsic and self.args.f_norm:
                            # full uncertainty
                            y_new = self.y_extrinsic_ma[dset][i_e] * self.f_inv[dset]
                            dy = (y_new - y_ma) / y_ma
                            self.ax_resid_shift.errorbar(E, dy.mean, yerr=dy.sdev,
                                                   marker='s', mfc='None', c=clr, alpha=0.5)
                            '''
                            # only extrinsic uncertainty
                            y_new = self.y_extrinsic_ma[dset][i_e] * self.f_inv[dset].mean
                            dy = (y_new - y_ma) / y_ma
                            self.ax_resid.errorbar(E, dy.mean, yerr=dy.sdev,
                                                   marker='None', c='k', alpha=0.5)
                            # only original stoch uncertainty
                            y_new = self.y[dset][i_e] * self.f_inv[dset]
                            dy = (y_new - y_ma) / y_ma
                            self.ax_resid.errorbar(E, dy.mean, yerr=dy.sdev,
                                                   marker='o', c=clr, mfc='None')
                            '''

                    self.ax_resid.set_xlim(plot_params['x_lim'])
                    self.ax_resid.set_ylim(-0.65, 1.05)
                    self.ax_resid.set_xscale(plot_params['x_scale'])
                    self.ax_resid.set_xlabel(plot_params['x_label'], fontsize=12)
                    self.ax_resid.set_ylabel(r'$\Delta S/S$', fontsize=12)

                    self.ax_resid_shift.set_xlim(plot_params['x_lim'])
                    self.ax_resid_shift.set_ylim(-0.65, 1.05)
                    self.ax_resid_shift.set_xscale(plot_params['x_scale'])
                    self.ax_resid_shift.set_xlabel(plot_params['x_label'], fontsize=12)


    def report_fits(self):
        for model in self.fit_results:
            print('-----------------------------------------------------------')
            print(model)
            print('-----------------------------------------------------------')
            print(self.fit_results[model])
