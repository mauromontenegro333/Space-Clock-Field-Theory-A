"""Compile the standalone integrated manuscript with pdfLaTeX three times."""
import argparse
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'build')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = ROOT / 'manuscript' / 'SCFT_Integrated.tex'
    shutil.copy2(source, out / source.name)
    for _ in range(3):
        subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',source.name],cwd=out,check=True)
    print(out / 'SCFT_Integrated.pdf')

if __name__ == '__main__':
    main()
