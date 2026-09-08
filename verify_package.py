"""Verify package hashes, numerical schemas, scan keys, and exact-check reports."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    manifest = ROOT / 'MANIFEST_SHA256.txt'
    checked = 0
    for line in manifest.read_text().splitlines():
        expected, relative = line.split('  ', 1)
        path = ROOT / relative
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f'hash mismatch: {relative}')
        checked += 1

    reference = np.loadtxt(ROOT / 'data/reference/scft_background.csv', delimiter=',')
    refinement = np.loadtxt(ROOT / 'data/refinement/scft_background.csv', delimiter=',')
    sourcefree = np.loadtxt(ROOT / 'data/sourcefree/scft_sourcefree.csv', delimiter=',')
    assert reference.shape == (20001, 24)
    assert refinement.shape == (40001, 24)
    assert sourcefree.shape == (4001, 24)

    with np.load(ROOT / 'data/reference/scft_modes.npz', allow_pickle=False) as modes:
        assert set(modes.files) == {'0.1', '1.0', '10.0', '100.0', '1000.0'}
        assert all(modes[key].shape == (6001, 7) for key in modes.files)
    with np.load(ROOT / 'data/scan/scft_modes.npz', allow_pickle=False) as modes:
        assert len(modes.files) == 42
        assert all(modes[key].shape == (6001, 7) for key in modes.files)
    with np.load(ROOT / 'data/refinement/scft_modes.npz', allow_pickle=False) as modes:
        assert set(modes.files) == {'10.0', '362.0291470135738'}
        assert all(modes[key].shape == (6001, 7) for key in modes.files)

    completed = json.loads((ROOT / 'checks/scft_completion_checks.json').read_text())
    rank1 = json.loads((ROOT / 'checks/scft_rank1_checks.json').read_text())
    assert completed['all_passed'] and completed['passed_count'] == completed['total_count'] == 14
    assert rank1['all_passed'] and rank1['passed_count'] == rank1['total_count'] == 54
    assert completed['floating_point_used'] is False
    assert rank1['floating_point_used'] is False
    print(json.dumps({'manifest_files_checked': checked,
                      'reference_background_shape': list(reference.shape),
                      'scan_modes': 42,
                      'completed_exact_checks': 14,
                      'rank1_exact_checks': 54,
                      'all_passed': True}, indent=2))

if __name__ == '__main__':
    main()
