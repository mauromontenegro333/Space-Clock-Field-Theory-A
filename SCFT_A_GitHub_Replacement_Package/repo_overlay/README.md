# Space-Clock Field Theory A — reproducibility and theorem-support package

This repository contains the reproducibility material for the revised SCFT-A classical construction together with the algebraic/support package for the local nonlinear vacuum Cauchy theorem.

The manuscript and technical-companion TeX/PDF files are intentionally **not distributed in this repository**. This repository is for scripts, exact certificates, numerical data, research ledgers, and theorem-support artifacts only.

## Revised classical benchmark

The exact revised coefficient point used by the current numerical package is

- `kappa_m / H_d^2 = 10^30`
- `delta = 1 / 10^140`
- `sigma = 1 / 10^6`
- `a_0 / H_d = 1`
- `s_N = 10`
- `Q_d = 103 / 100`
- `b_d / H_d = 80 / 103`
- `kappa_d / H_d^2 = -(80 / 103)^2 / 10`
- `nu_A = 30`
- `lambda_D = 1`

The coefficient certificate records 33 passed sufficient-bound checks. The exact algebra package records 16 passed algebra/scale checks. The numerical package recomputes the source-free and visible-loaded backgrounds, the endpoint spectrum, a contracting branch, and 42 finite-time linear responses.

For the stored response calculation, the maximum sampled unit response is `1.0011357323900973`, the maximum channel-scaled refinement difference is `1.6975908959798724e-11`, the maximum slip residual is `7.860465750519907e-19`, and the minimum canonical denominator is `1.792669897832447`.

## Local nonlinear vacuum Cauchy result

The directory `nonlinear_wellposedness/` contains the machine-readable scope certificate, exact symbol checks, exact benchmark values, operator-order ledger, literature/search ledgers, and the symbolic verifier associated with the newer local nonlinear well-posedness result.

On the stated open strongly Dirac-regular component of the compact boundaryless three-torus, for integer `s >= 18`, the certified theorem scope is local Hadamard well-posedness in `H^(s+3) x H^s`, including local existence, uniqueness, continuous dependence, persistence of higher regularity, constraint propagation, and a continuation criterion.

The theorem scope does **not** claim global existence, long-time nonlinear stability, a bounded-domain initial-boundary-value theorem, unconditional evolution for unspecified `S_vis`, or quantum/strong-coupling control. The symbolic verifier certifies exact algebra and artifact integrity; it does not machine-prove the functional-analytic theorem.

## Reproduce the revised classical calculations

Use Python 3.12 and the pinned dependencies:

```bash
python3 -m pip install -r requirements.txt
python3 scripts/certify_coefficients.py
python3 scripts/exact_algebra.py
python3 scripts/reproduce_numerics.py
```

The exact coefficient and algebraic scripts use rational/symbolic or rigorously enclosed quantities for their stated checks. CSV decimal columns and ODE outputs are numerical approximations.

## Reproduce the nonlinear algebraic certificate

```bash
cd nonlinear_wellposedness
python3 -m pip install -r requirements.txt
python3 scripts/verify_wellposedness.py
```

The proof report itself is intentionally not bundled here.

## Repository layout

- `scripts/` — revised coefficient certification, exact algebra, and numerical reproduction.
- `data/` — revised benchmark certificates, backgrounds, endpoint/contracting data, response summaries, representative time series, and the compressed full response array.
- `nonlinear_wellposedness/` — support artifacts for the local nonlinear vacuum Cauchy theorem.
- `REVISION_REPORT.md` — current status of the eight major theoretical/reproducibility issues.
- `SHA256SUMS.txt` — hashes of the distributed repository files.
- `CITATION.cff` and `ZENODO_METADATA.md` — current package metadata.

## Scope

This is a classical theory/reproducibility repository. It does not provide a CMB likelihood, primordial Boltzmann spectrum, nonlinear galaxy matching, a complete matter strong-coupling domain, radiative protection, or a quantum completion.
