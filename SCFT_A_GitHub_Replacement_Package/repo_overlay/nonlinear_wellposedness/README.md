# SCFT-A local nonlinear well-posedness — support artifacts

This directory contains the algebraic, numerical-free, and research-ledger support artifacts associated with the newer local nonlinear vacuum Cauchy theorem for the completed SCFT-A action.

The theorem report TeX/PDF is intentionally **not bundled in this repository**. The machine-readable certificate records the theorem scope; the symbolic verifier checks exact finite-dimensional identities and artifact integrity only.

## Certified theorem scope

On an open strongly Dirac-regular component of a compact boundaryless three-torus, for integer `s >= 18`, the theorem scope is local Hadamard well-posedness in `H^(s+3) x H^s`, including:

- local existence;
- uniqueness;
- continuous dependence;
- persistence of higher regularity;
- constraint propagation;
- a continuation criterion.

Excluded claims are global existence, long-time nonlinear stability, bounded-domain initial-boundary-value well-posedness, unconditional evolution for unspecified `S_vis`, and quantum/strong-coupling control.

## Reproduce the exact algebraic certificate

```bash
python3 -m pip install -r requirements.txt
python3 scripts/verify_wellposedness.py
```

The verifier checks the exact principal symmetrizer, frequency multiplicities, positive Hessian spectra, the second-class constraint inverse (including a noncommuting block check), the lapse Schur sign, and the exact rolling benchmark identity/signs.

## Files

- `scripts/verify_wellposedness.py` — exact symbolic verifier and manifest generator.
- `data/proof_certificate.json` — machine-readable theorem scope and PASS results.
- `data/exact_symbol_checks.json` — exact matrices, spectra, and zero residuals.
- `data/exact_benchmark.csv` — exact unrounded representative-point values.
- `data/operator_order_ledger.csv` — operator-order accounting.
- `data/literature_sources.csv` — primary-source applicability ledger.
- `data/search_queries.csv` — reproducible search record.
- `data/sha256_manifest.txt` — hashes of this distributed support package.

The symbolic verifier does not replace expert review of the functional-analytic proof.
