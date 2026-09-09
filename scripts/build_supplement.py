#!/usr/bin/env python3
from pathlib import Path
import re, shutil, subprocess

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "manuscript" / "SCFT_A_Technical_Companion.tex"
OUT = ROOT / "manuscript" / "SCFT_A_Technical_Companion.pdf"
BUILD = ROOT / "build"
BUILD.mkdir(exist_ok=True)
engine = shutil.which("pdflatex")
if engine is None:
    raise SystemExit("pdflatex is required (TeX Live or MiKTeX).")
for i in range(3):
    result = subprocess.run([
        engine, "-interaction=nonstopmode", "-halt-on-error", "-file-line-error",
        "-output-directory", str(BUILD), str(SRC)
    ], cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (BUILD / f"supplement_pass_{i+1}.log").write_text(result.stdout)
    if result.returncode:
        raise SystemExit(result.stdout[-7000:])
log = (BUILD / "SCFT_A_Technical_Companion.log").read_text(errors="replace")
bad = [p for p in [
    r"There were undefined references", r"There were undefined citations",
    r"Overfull \\hbox", r"Overfull \\vbox", r"! LaTeX Error"
] if re.search(p, log)]
if bad:
    raise SystemExit("Build audit failed: " + ", ".join(bad))
shutil.copy2(BUILD / "SCFT_A_Technical_Companion.pdf", OUT)
print(OUT)
