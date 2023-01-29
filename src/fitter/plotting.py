import matplotlib.pyplot as plt


def plot_data(x,y, plot_struct):
    plt.figure(figsize=(7,4))
    ax = plt.axes([0.12,0.12,0.86,0.86])

    for dset in x:
        y_m = [k.mean for k in y[dset]]
        y_s = [k.sdev for k in y[dset]]
        ax.errorbar(x[dset], y_m, yerr=y_s,
                    marker='o', mfc='None', linestyle='None',
                    label=dset.replace('_','\_'))
    ax.set_xlim(plot_struct['x_lim'])
    ax.set_ylim(plot_struct['y_lim'])
    ax.set_xscale(plot_struct['x_scale'])
    ax.set_yscale(plot_struct['y_scale'])
    ax.set_xlabel(plot_struct['x_label'], fontsize=16)
    ax.set_ylabel(plot_struct['y_label'], fontsize=16)
    ax.legend(loc=4)

    return ax
