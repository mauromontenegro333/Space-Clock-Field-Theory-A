"""Run both exact checkers in a new directory and save their reports."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'exact_check_rerun')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    subprocess.run([sys.executable, str(ROOT / 'scripts/scft_completion_checks.py')], cwd=out, check=True)
    subprocess.run([sys.executable, str(ROOT / 'scripts/scft_rank1_checks.py'),
                    '--output', str(out / 'scft_rank1_checks.json')], cwd=out, check=True)
    reports = [json.loads((out / name).read_text()) for name in
               ['scft_completion_checks.json', 'scft_rank1_checks.json']]
    assert [(r['passed_count'], r['total_count'], r['all_passed']) for r in reports] == [
        (14, 14, True), (54, 54, True)]
    print(out)

if __name__ == '__main__':
    main()
