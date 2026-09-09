# Space-Clock Field Theory A — Technical Supplement & Reproducibility

Publication support repository for **SCFT-A: Technical Companion and Reproducibility Certificate** by Mauro Alfonso Montenegro.

This repository intentionally contains **only the technical supplement and its reproducibility materials**. The long-form main manuscript is not included here in either PDF or LaTeX form.

## Publication file
- `manuscript/SCFT_A_Technical_Companion.pdf` — final 31-page technical supplement.
- `manuscript/SCFT_A_Technical_Companion.tex` — exact integrated LaTeX source.

## Reproducibility
- `scripts/verify_wellposedness.py` — exact algebraic/structural verification for the local nonlinear vacuum Cauchy theorem.
- `data/proof_certificate.json` — theorem scope and PASS ledger.
- `data/exact_symbol_checks.json` — exact matrices, spectra, Schur sign, and zero residuals.
- `data/exact_benchmark.csv` — exact unrounded benchmark values.
- `data/operator_order_ledger.csv` — operator-order accounting.

The local theorem concerns the vacuum gravity-clock sector on an open strongly Dirac-regular component of a compact boundaryless three-torus. It does not claim global nonlinear stability, bounded-domain well-posedness, unrestricted visible-matter evolution, or quantum/strong-coupling control. Machine checks certify algebra and artifact integrity; the functional-analytic proof remains the technical supplement itself.

## Build
```bash
python3 -m pip install -r requirements.txt
python3 scripts/verify_wellposedness.py
```
To build the supplement, run `pdflatex` three times on `manuscript/SCFT_A_Technical_Companion.tex`.

## License
CC0 1.0 Universal. See `LICENSE`.
