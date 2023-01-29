import numpy as np
import gvar as gv

import fitter.load_data as ld
import fitter.plotting as plot

class FitManager():

    def __init__(self,args):
        self.x, self.y, self.plot_struct = ld.load_data(args.d_file, args.d_path, args.d_sets)

    def plot_data(self):
        self.ax = plot.plot_data(self.x,self.y, self.plot_struct)
