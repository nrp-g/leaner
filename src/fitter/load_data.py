import tables as h5
import numpy as np
import gvar as gv

def load_data(d_file, d_path, dsets=[]):

    y = {}
    x = {}
    plot_struct = {}

    with h5.open_file(d_file,'r') as f5:
        for attr in ['x_label', 'y_label', 'x_scale', 'y_scale', 'x_lim', 'y_lim']:
            if attr in f5.get_node('/'+d_path)._v_attrs:
                plot_struct[attr] = f5.get_node('/'+d_path)._v_attrs[attr]
        # if dsets is non empty list
        if dsets:
            data_sets = dsets
        else:
            data_sets = []
            for g in f5.get_node('/'+d_path):
                data_sets.append(g._v_name)
        for dset in data_sets:
            if 'dS' in f5.get_node('/'+d_path+'/'+dset):
                x[dset] = f5.get_node('/'+d_path+'/'+dset+'/E').read()
                tmp_y   = f5.get_node('/'+d_path+'/'+dset+'/S').read()
                tmp_dy  = f5.get_node('/'+d_path+'/'+dset+'/dS').read()
                y[dset] = gv.gvar(tmp_y, tmp_dy)

    return x,y, plot_struct
