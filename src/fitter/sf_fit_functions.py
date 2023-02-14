import numpy as np
import gvar as gv

class SFFunctions():

    def __init__(self, model, f_norm=True, offset=True, pheno_file='data/Spd_newrun_interspline.dat'):
        self.model  = model.split('_')[0]
        try:
            self.order = int(model.split('_')[-1])
        except:
            self.order = 0
        self.f_norm = f_norm
        self.offset = offset

        if 'pheno' in model:
            from scipy.interpolate import interp1d
            pheno_spline = open(pheno_file).readlines()
            e_mev = []
            S     = []
            for l in pheno_spline:
                e_mev.append(float(l.split()[0])/1000)
                S.append(float(l.split()[1])*1e-6)
            e_mev = np.array(e_mev)
            S     = np.array(S)
            self.S_pheno = interp1d(e_mev,S, kind='cubic')


    def fit_func(self, x, p):
        # use the model name to determine the fit function to use
        return getattr(self, self.model)(x, p)

    def prune_priors(self, priors, dsets):
        # prune the priors so we only use ones that are used
        priors_pruned = {}
        p0            = {}
        for p in priors:
            if 'log' not in p:
                if 'pheno' not in self.model and 'S' in p:
                    if int(p.split('_')[-1]) <= self.order:
                        priors_pruned[p] = gv.gvar(priors[p].mean, priors[p].sdev)
                        p0[p]            = priors[p].mean
                else:
                    if (self.offset and p in ['a', 'b']) or (p == 'a'):
                        priors_pruned[p] = gv.gvar(priors[p].mean, priors[p].sdev)
                        p0[p]            = priors[p].mean

            elif 'log' in p and self.f_norm:
                if any(k in p for k in dsets):
                    priors_pruned[p] = gv.gvar(priors[p].mean, priors[p].sdev)
                    p0[p]            = priors[p].mean

        return priors_pruned, p0

    def poly(self, x, p):
        result = {}
        for dset in x:
            result[dset] = p['S_0']
            for n in range(1, self.order+1):
                result[dset] += p['S_%d' %n] * x[dset]**n
            if self.f_norm:
                result[dset] = result[dset] * p['f_'+dset]

        return result

    def pheno(self, x, p):
        result = {}
        for dset in x:
            result[dset] = p['a'] * self.S_pheno(x[dset])
            if self.offset:
                result[dset] += p['b']
            if self.f_norm:
                result[dset] = result[dset] * p['f_'+dset]

        return result
