import gvar as gv

class SFFunctions():

    def __init__(self, model, f_norm=True):
        self.model  = model.split('_')[0]
        self.order  = int(model.split('_')[-1])
        self.f_norm = f_norm

    def fit_func(self, x, p):
        # use the model name to determine the fit function to use
        return getattr(self, self.model)(x, p)

    def prune_priors(self, priors, dsets):
        # prune the priors so we only use ones that are used
        priors_pruned = {}
        p0            = {}
        for p in priors:
            if 'log' not in p:
                if int(p.split('_')[-1]) <= self.order:
                    priors_pruned[p] = gv.gvar(priors[p].mean, priors[p].sdev)
                    p0[p]            = priors[p].mean
            elif 'log' in p and self.f_norm:
                if any(k in p for k in dsets):
                    priors_pruned[p] = gv.gvar(priors[p].mean, priors[p].sdev)
                    p0[p]            = priors[p].mean

        return priors_pruned, p0

    def poly(self, x, p):
        y = {}
        for dset in x:
            y[dset] = p['S_0']
            for n in range(1, self.order+1):
                y[dset] += p['S_%d' %n] * x[dset]**n
            if self.f_norm:
                y[dset] = y[dset] * p['f_'+dset]
        return y
