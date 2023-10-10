import matplotlib.pyplot as plt


def plot_data(x,y, plot_params, marker='o', color=True, label=True, error=True, alpha=0.4):
    if plot_params['residuals']:
        plt.figure(figsize=(7,5))
        ax = plt.axes([0.13,0.3,0.85,0.68])
        ax_resid = plt.axes([0.13,0.13,0.85,0.17])
    else:
        plt.figure(figsize=(7,4))
        ax = plt.axes([0.13,0.13,0.85,0.85])

    for dset in x:
        y_m = [k.mean for k in y[dset]]
        y_s = [k.sdev for k in y[dset]]
        if color:
            clr = plot_params['colors'][dset]
        else:
            clr = 'k'
        if label:
            lbl = dset.replace('_','\_')
        else:
            lbl = ''
        if error:
            ax.errorbar(x[dset], y_m, yerr=y_s, alpha=alpha,
                        marker=marker, c=clr, mfc='None', linestyle='None', label=lbl)
        else:
            ax.plot(x[dset], y_m, alpha=alpha,
                    marker=marker, c=clr, mfc='None', linestyle='None', label=lbl)

    ax.set_xlim(plot_params['x_lim'])
    ax.set_ylim(plot_params['y_lim'])
    ax.set_xscale(plot_params['x_scale'])
    ax.set_yscale(plot_params['y_scale'])
    if not plot_params['residuals']:
        ax.set_xlabel(plot_params['x_label'], fontsize=16)
    else:
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.tick_params(direction='inout', which='both')

    ax.set_ylabel(plot_params['y_label'], fontsize=16)
    if label:
        if 'legend_loc' in plot_params:
            ax.legend(loc=plot_params['legend_loc'])
        else:
            ax.legend()
    ax.text(plot_params['text_loc'][0], plot_params['text_loc'][1],
            plot_params['text'], transform=ax.transAxes,fontsize=20)

    if plot_params['residuals']:
        return ax, ax_resid
    else:
        return ax
