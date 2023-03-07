import sys, os
import argparse
#try:
import fitter.fit_manager as FM
import fitter.load_data as ld
#except:
#    print('running without installation')
#    print('THIS DOES NOT WORK YET')
#    import src.fitter.fit_manager as FM

def main():
    parser = argparse.ArgumentParser(
        description='Perform Bayes analysis of nuclear reaction data S-factors')
    parser.add_argument('--SF_reaction',    type=str, default='S_12',
                        help=               'specify SF reaction [%(default)s]'  )
  
    parser.add_argument('--fit_params',     default='input/fit_params.py',
                        help=               'specify user input file [%(default)s]')

    parser.add_argument('--d_file',         default='data/solar_fusion_reacions.h5',
                        help=               'specify data file [%(default)s]')
    parser.add_argument('--d_sets',         nargs='+',
                        help=               'user can specify d_set list to override input file')

    parser.add_argument('--models',         nargs='+', help='overide models in input file with this list')
    parser.add_argument('--f_norm',         default=True, action='store_false',
                        help=               'S -> f S where f is an normalization'
                                            +'factor to be determined [%(default)s]')
    parser.add_argument('--extrinsic',      nargs='+', default=['rel'],
                        help=               'list of extrinsic statistical uncertainty models'
                                            +"in analysis, options are rel, abs and '' [%(default)s]")
    parser.add_argument('--pheno_file',     default='data/Marcucci_etal_PRL116_2016_1510.07877.dat',
                        help=               'what .dat file to use for pheno function input [%(default)s]')

    parser.add_argument('--run_analysis',   default=True, action='store_false',
                        help=               'run Bayes Model Analysis? [%(default)s]')
    parser.add_argument('--redo_fits',      default=False, action='store_true',
                        help=               'redo fits even if saved? [%(default)s]')
    parser.add_argument('--report_fits',    default=False, action='store_true',
                        help=               'print results from each model [%(default)s]')
    parser.add_argument('--report_models',  default=True, action='store_true',
                        help=               'report model weights? [%(default)s]')

    parser.add_argument('--prior_width',    nargs='+', type=str,
                        help=               'list priors to simultaneously perform a width study')
    parser.add_argument('--prior_range',    nargs='+', type=float,
                        help=               'min, max, ds: to be used as np.arange(min,max+dm,dm) '
                                            +'to be used for prior_width study')

    parser.add_argument('--show_plot',      default=True, action='store_false',
                        help=               'show plots? [%(default)s]')

    parser.add_argument('--interact',       default=False, action='store_true',
                        help=               'open IPython instance after to interact with results? [%(default)s]')

    args = parser.parse_args()
    print(args)

    # check for problems
    # turn this into a function that checks everything
    if (args.prior_width and not args.prior_range) or (not args.prior_width and args.prior_range):
        sys.exit('you must specify BOTH prior_width and prior_range to perform prior_width study')

    # populate the args with essential info
    #args.d_sets = []
    args.d_path = args.SF_reaction

    if args.show_plot:
        import matplotlib.pyplot as plt
        plt.ion()

    # start the fit manager
    sf_fit = FM.FitManager(args)

    # perform prior_width study
    if args.prior_width:
        sf_fit.prior_width()

    if args.run_analysis:
        # perform fit over models
        sf_fit.fit_models()

        # report fits
        if args.report_fits:
            sf_fit.report_fits()

        if args.report_models:
            sf_fit.report_models()

        sf_fit.model_avg_S(0.0,   plot_hist=True)
        sf_fit.model_avg_S(0.091, plot_hist=True)

        # plot data
        if args.show_plot:
            sf_fit.plot_data()
            sf_fit.plot_fit()

            if not os.path.exists('figures'):
                os.makedirs('figures')
            plt.savefig('figures/S_model_avg.pdf', transparent=True)

            #sf_fit.plot_pdf()

    if args.show_plot:
        plt.ioff()
        plt.show()

    if args.interact:
        import IPython; IPython.embed()


if __name__ == "__main__":
    main()
