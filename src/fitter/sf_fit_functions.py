import gvar as gv

class SFFunctions():

    def __init__(self, model):
        self.model = model.split()[0]
        self.order = int(model.split()[-1])

    def fit_func(self, x, p):
        # use the model name to determine the fit function to use
        return getattr(self, self.model)(x, p)

    def prune_priors(self, priors):
        # prune the priors so we only use ones that are used
        priors_pruned = {}
        p0            = {}
        for p in priors:
            if int(p.split('_')[-1]) <= self.order:
                priors_pruned[p] = gv.gvar(priors[p].mean, priors[p].sdev)
                p0[p]            = priors[p].mean
        return priors_pruned, p0

    def poly(self, x, p):
        y = {}
        for dset in x:
            y[dset] = p['S_0']
            for n in range(1, self.order+1):
                y[dset] += p['S_%d' %n] * x[dset]**n
        return y
