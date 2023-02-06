import matplotlib.pyplot as plt


def plot_data(x,y, plot_params):
    plt.figure(figsize=(7,4))
    ax = plt.axes([0.13,0.13,0.85,0.85])

    for dset in x:
        y_m = [k.mean for k in y[dset]]
        y_s = [k.sdev for k in y[dset]]
        clr = plot_params['colors'][dset]
        ax.errorbar(x[dset], y_m, yerr=y_s,
                    marker='o', c=clr, mfc='None', linestyle='None',
                    label=dset.replace('_','\_'))
    ax.set_xlim(plot_params['x_lim'])
    ax.set_ylim(plot_params['y_lim'])
    ax.set_xscale(plot_params['x_scale'])
    ax.set_yscale(plot_params['y_scale'])
    ax.set_xlabel(plot_params['x_label'], fontsize=16)
    ax.set_ylabel(plot_params['y_label'], fontsize=16)
    if 'legend_loc' in plot_params:
        ax.legend(loc=plot_params['legend_loc'])
    else:
        ax.legend()
    ax.text(plot_params['text_loc'][0], plot_params['text_loc'][1],
            plot_params['text'], transform=ax.transAxes,fontsize=20)

    return ax
