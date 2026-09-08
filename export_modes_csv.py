"""Export NPZ modes to labeled CSV without reducing stored float precision."""
import argparse
from pathlib import Path
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    with np.load(args.input, allow_pickle=False) as modes:
        for key in modes.files:
            # Parse the stored key before using it in a filename.
            kval = float(key)
            name = 'k_' + repr(kval).replace('.', 'p') + '.csv'
            np.savetxt(args.output / name, modes[key], delimiter=',', fmt='%.18e',
                       header='N,zeta,Psi_N,delta_b_N,delta_r_N,P_over_H2,k_over_aH')
            print(key, name)


if __name__ == '__main__':
    main()
