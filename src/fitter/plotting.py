import matplotlib.pyplot as plt


def plot_data(x,y, args):
    plt.figure(figsize=(7,4))
    ax = plt.axes([0.12,0.12,0.86,0.86])

    for dset in x:
        y_m = [k.mean for k in y[dset]]
        y_s = [k.sdev for k in y[dset]]
        ax.errorbar(x[dset], y_m, yerr=y_s,
                    marker='o', mfc='None', linestyle='None',
                    label=dset.replace('_','\_'))
    ax.set_xlim(6.5e-4, 2.5)
    ax.set_ylim(8e-8,   2.5e-5)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'Energy (MeV)',fontsize=16)
    ax.set_ylabel(r'S-factor (MeV b)',fontsize=16)
    ax.legend(loc=4)

    return ax
