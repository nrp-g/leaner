import sys
import argparse
try:
    import fitter.fit_manager as FM
    import fitter.load_data as ld
except:
    print('running without installation')
    print('THIS DOES NOT WORK YET')
    import src.fitter.fit_manager as FM

def main():
    parser = argparse.ArgumentParser(
        description='Perform Bayes analysis of nuclear reaction data S-factors')
    parser.add_argument('--SF_reaction', type=str, default='S_12',
                        help=            'specify SF reaction [%(default)s]')

    parser.add_argument('--d_file',      default='data/solar_fusion_reacions.h5',
                        help=            'specify data file [%(default)s]')

    parser.add_argument('--show_plot',   default=True, action='store_false',
                        help=            'show plots? [%(default)s]')

    parser.add_argument('--interact',    default=False, action='store_true',
                        help=            'open IPython instance after to interact with results? [%(default)s]')

    args = parser.parse_args()
    print(args)

    # populate the args with essential info
    args.d_sets = []
    args.d_path = args.SF_reaction

    # start the fit manager
    sf_fit = FM.FitManager(args)

    if args.show_plot:
        import matplotlib.pyplot as plt
        plt.ion()
        sf_fit.plot_data()

    if args.interact:
        import IPython; IPython.embed()

    if args.show_plot:
        plt.ioff()
        plt.show()

if __name__ == "__main__":
    main()
