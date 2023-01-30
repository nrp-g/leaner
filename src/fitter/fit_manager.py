import importlib
import sys
import os
import numpy as np
import gvar as gv

import fitter.load_data as ld
import fitter.plotting as plot

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
