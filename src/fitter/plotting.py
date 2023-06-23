import matplotlib.pyplot as plt


def plot_data(x,y, plot_params, marker='o', color=True, label=True, error=True, alpha=0.4):
    if plot_params['residuals']:
        plt.figure(figsize=(6.5,4))
        ax = plt.axes([0.11,0.3,0.44,0.68])
        ax_resid = plt.axes([0.11,0.11,0.44,0.17])
        ax_shift = plt.axes([0.55,0.3,0.44,0.68])
        ax_resid_shift = plt.axes([0.55,0.11,0.44,0.17])
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
        ax.set_xlabel(plot_params['x_label'], fontsize=12)
    else:
        ax_shift.set_xlim(plot_params['x_lim'])
        ax_shift.set_ylim(plot_params['y_lim'])
        ax_shift.set_xscale(plot_params['x_scale'])
        ax_shift.set_yscale(plot_params['y_scale'])

        ax.xaxis.set_tick_params(labelbottom=False)
        ax.tick_params(direction='inout', which='both', right=True)

        ax_shift.xaxis.set_tick_params(labelbottom=False)
        ax_shift.yaxis.set_tick_params(labelleft=False)
        ax_resid_shift.yaxis.set_tick_params(labelleft=False)
        ax_shift.tick_params(direction='inout', which='both', right=True)

    ax.set_ylabel(plot_params['y_label'], fontsize=12)
    if label:
        if 'legend_loc' in plot_params:
            ax.legend(loc=plot_params['legend_loc'])
        else:
            ax.legend()
    ax.text(plot_params['text_loc'][0], plot_params['text_loc'][1],
            plot_params['text'], transform=ax.transAxes,fontsize=12)
    ax_shift.text(0.1, 0.8, r'$y_{D_i}^{\rm ext} / f_D$', transform=ax_shift.transAxes,fontsize=12)

    if plot_params['residuals']:
        return ax, ax_resid, ax_shift, ax_resid_shift
    else:
        return ax
