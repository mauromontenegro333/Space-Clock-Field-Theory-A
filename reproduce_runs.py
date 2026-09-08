"""Rerun the base-action S_0 numerical profiles represented by the delivered data files.
These runs do not compute perturbations of the completed spatial action."""
import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import scft_audit as audit


def run_profile(profile, destination):
    destination.mkdir(parents=True, exist_ok=False)
    previous = Path.cwd()
    try:
        os.chdir(destination)
        if profile == 'reference':
            audit.run(ks=[.1, 1., 10., 100., 1000.])
        elif profile == 'scan':
            audit.run(ks=np.geomspace(.03, 1000, 42))
        elif profile == 'refinement':
            audit.run(npoints=40001, ks=[10., 362.0291470135738],
                      rtol=2e-10, atol=2e-14, max_step=.005)
        elif profile == 'sourcefree':
            ns = np.linspace(np.log(1e-10), np.log(900), 4001)
            bg = np.array([audit.background(n, loaded=False) for n in ns])
            np.savetxt('scft_sourcefree.csv', np.column_stack([ns, bg]), delimiter=',')
    finally:
        os.chdir(previous)


def compare_refinement(destination):
    differences = {}
    with np.load(destination / 'refinement/scft_modes.npz', allow_pickle=False) as fine:
        for key, baseline in [('10.0', 'reference'),
                              ('362.0291470135738', 'scan')]:
            with np.load(destination / baseline / 'scft_modes.npz', allow_pickle=False) as base:
                x, y = fine[key], base[key]
                if not np.array_equal(x[:, 0], y[:, 0]):
                    raise ValueError('Output time grids differ; no implicit interpolation applied.')
                differences[key] = {
                    name: float(np.max(np.abs(x[:, i] - y[:, i]) /
                                       np.maximum(1., np.abs(x[:, i]))))
                    for i, name in enumerate(['zeta', 'Psi', 'delta_b', 'delta_r'], start=1)
                }
    (destination / 'refinement_differences.json').write_text(
        json.dumps(differences, indent=2) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('profile', choices=['reference', 'scan', 'refinement', 'sourcefree', 'all'])
    parser.add_argument('--output', type=Path, default=ROOT / 'rerun')
    args = parser.parse_args()
    destination = args.output.resolve()
    profiles = ['reference', 'scan', 'refinement', 'sourcefree'] if args.profile == 'all' else [args.profile]
    for profile in profiles:
        print(f'Running {profile} in {destination / profile}', flush=True)
        run_profile(profile, destination / profile)
    if args.profile == 'all':
        compare_refinement(destination)


if __name__ == '__main__':
    main()
