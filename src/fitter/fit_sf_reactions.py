import sys
import argparse
try:
    import fitter.fit_manager as FM
except:
    print('running without installation')
    print('THIS DOES NOT WORK YET')
    import src.fitter.fit_manager as FM

def main():
    parser = argparse.ArgumentParser(
        description='Perform Bayes analysis of nuclear reaction data S-factors')
    parser.add_argument('--SF_reaction', type=str, default='S_12',
                            help=            'specify SF reaction [%(default)s]')


    args = parser.parse_args()
    print(args)

if __name__ == "__main__":
    main()
