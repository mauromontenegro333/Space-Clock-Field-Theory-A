# SCFT current manuscript data and code

This archive accompanies the 8 September 2026 integrated manuscript:

**Space-Clock Field Theory: A Single-Metric Theory of Dark-Sector Unification**  
Mauro Alfonso Montenegro

The package contains only files still used by, quoted in, or required to reproduce retained calculations in the current 92-page manuscript. Exact duplicates, older manuscript assembly scripts, obsolete action repairs, external copies of figures already embedded in the LaTeX source, and build logs have been excluded.

## Scientific scope

The package separates two scopes that must not be conflated:

1. `data/` and `scripts/scft_audit.py` support the retained **base-action `S_0`** homogeneous and finite-wavelength results. The current manuscript labels those figures, transfer functions, and mode scans as `S_0` calculations. They are not finite-momentum predictions of the completed spatial action.
2. `scripts/scft_completion_checks.py`, `scripts/scft_rank1_checks.py`, and `checks/` support exact algebraic identities for the **completed local transition**. They establish the reported quadratic and auxiliary algebra. They do not prove the loss-free nonlinear energy estimate or global constraint compatibility. Full inhomogeneous nonlinear Rank-1 continuation remains conditional.

`FILE_STATUS.csv` gives the scope and role of each file group. `DATA_DICTIONARY.md` defines all stored columns and summary fields. `PROVENANCE.json` records the environment and proof-status boundary.

## Contents

- `manuscript/SCFT_Integrated.tex` — standalone LaTeX source with vector plot coordinates and source-code blocks embedded.
- `manuscript/SCFT_Integrated.pdf` — verified 92-page PDF.
- `scripts/scft_audit.py` — primary `S_0` numerical implementation.
- `scripts/reproduce_runs.py` — portable launcher for the four stored numerical profiles.
- `scripts/export_modes_csv.py` — precision-preserving NPZ exporter.
- `scripts/scft_completion_checks.py` — 14 exact completed-action checks.
- `scripts/scft_rank1_checks.py` — 54 exact local-transition checks.
- `scripts/historical_base_check/scft_local_vacuum_check.py` — appendix-level historical base diagnostic, explicitly not a completed-action check.
- `data/reference/` — 20,001-point matter-loaded background, five reference modes, and summary.
- `data/sourcefree/` — 4,001-point source-free background.
- `data/scan/` — complete 42-mode scan arrays and summary.
- `data/refinement/` — 40,001-point background, two tighter-tolerance modes, and convergence differences.
- `data/validation/` — refined extrema, regulator values, and one saved historical symbolic-result record.
- `checks/` — saved exact-check reports.

The reference and scan runs use the same stored 20,001-point background. Its exact duplicate was omitted.

## Python environment

The recorded environment is Python 3.12.13 with NumPy 2.3.5, SciPy 1.17.0, and SymPy 1.14.0.

```bash
python -m pip install -r requirements.txt
```

## Run the exact algebra checks

From the extracted package root:

```bash
python scripts/run_exact_checks.py --output exact_check_rerun
```

Expected results are `14/14` and `54/54`. These use symbolic/exact arithmetic and report `floating_point_used: false`.

## Reproduce the base-action numerical profiles

Each command writes into a new output directory and refuses to overwrite an existing profile directory.

```bash
python scripts/reproduce_runs.py reference --output rerun
python scripts/reproduce_runs.py scan --output rerun
python scripts/reproduce_runs.py refinement --output rerun
python scripts/reproduce_runs.py sourcefree --output rerun
```

To run all four profiles and reconstruct the refinement-difference report:

```bash
python scripts/reproduce_runs.py all --output complete_rerun
```

The scan and refinement integrations can take substantial computation time. Platform and library differences may affect the last binary64 digits and solver work counts.

## Export NPZ mode data to CSV

```bash
python scripts/export_modes_csv.py data/scan/scft_modes.npz exported_scan_csv
```

The exporter prints 18 digits after the decimal point so stored binary64 values round-trip through text.

## Compile the manuscript

Install a TeX distribution containing REVTeX 4.2, PGFPlots, Latin Modern, and the packages declared in the source. Then run:

```bash
python scripts/build_manuscript.py --output manuscript_build
```

The script runs pdfLaTeX three times. The LaTeX source also writes embedded Python sources through `filecontents*` when compiled.

## Integrity check

```bash
python scripts/verify_package.py
```

This validates `MANIFEST_SHA256.txt`, stored schemas, scan keys, and exact-check report counts.

## Interpretation limits

The mode integrations use minimally coupled barotropic baryon, photon, and neutrino perfect fluids. They are not a recombination calculation, collisionless-neutrino hierarchy, CMB or matter-power likelihood analysis, or completed-action perturbation scan. Numerical outputs are binary64 approximations; printed decimal digits do not make them exact mathematical constants.

No license has been assigned. Select an appropriate license in Zenodo before publishing.
