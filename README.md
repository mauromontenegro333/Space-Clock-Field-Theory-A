# SCFT-A intrinsic spatial curvature completion

This repository contains the current reproducibility package for the intrinsic-spatial-curvature completion of Space-Clock Field Theory (SCFT-A).

The pre-completion base-action `S_0` repository state has been preserved on the branch `legacy-s0-2026-09-08`. Files on `main` are the current curvature-completed calculation package and should not be mixed with the older `S_0` finite-wavelength outputs.

## Exact action addition and positive result

The completed action adds

`- eta_R (^(3)R)^2`

inside the `M_Pl^2/2` action bracket, with

`eta_R = 1/(2 kappa_m)`.

For `kappa_m/H_d^2 = 120`,

`eta_R H_d^2 = 1/240`.

The new scalar restoring contribution is `+8 eta_R p^4`.

The exact endpoint proof covers every `p^2 > 0` and the entire interval `0 < delta < 10^-6`. It establishes `D_d > 0`, `K_d > 0`, `V_d > 0`, and `V_d(120)/H_d^2 > 460`. All 16 implemented exact checks pass. The source-free local-transition and exact Minkowski results remain positive. This endpoint statement does not prove global matter-loaded nonlinear stability.

## Reproduce

Python 3 with NumPy, SciPy, and SymPy is required.

```bash
python3 -m pip install -r requirements.txt
python3 scft_curvature_completion.py --output reproduced_results
```

The script reruns the 16 exact checks, reconstructs the loaded background, evolves all 42 modes at two resolutions, performs three independent RK45 comparisons, and writes the current numerical products. Results may differ in the last printed binary64 digits across platforms.

No observational input or online data is required.

## Numerical conventions

The dimensionless numerical equations use `H_d = M_Pl = 1`. This does **not** assert the physical ratio `H_d/M_Pl = 1`.

The numerical response interval is `1 <= a <= 10`, with 1001 common output times. The loaded background contains 8001 nodes.

The inherited decimal value of `delta` in `coupling_input.json` is an uncertified numerical evaluation printed in the supplied SCFT-A source. It is not a newly certified global minimization, and no decimal approximation enters the exact endpoint sign proof.

## Files

- `scft_curvature_completion.py` — complete generator for the current exact checks and numerical response package.
- `exact_checks.json` — 16 executed exact checks and symbolic endpoint formulas.
- `numerical_summary.json` — convergence, residual, response-maximum, and independent-integrator summary.
- `coupling_input.json` — inherited decimal input and its provenance.
- `loaded_background.csv` — 8001-node reconstructed matter-loaded background.
- `endpoint_spectrum.csv` — corrected endpoint kinetic/restoring spectrum used for plotting.
- `response_summary.csv` — summary of all 42 modes.
- `response_00.csv`, `response_24.csv`, `response_30.csv`, `response_35.csv`, `response_38.csv`, `response_41.csv` — representative fine-response time series.
- `all_responses.npz` — all 42 fine responses; `values.shape == (42, 1001, 15)`.
- `inherited_background_figure.csv` — unchanged broad-range homogeneous figure coordinates extracted from the supplied TeX; not recomputed by the current response script.
- `DATA_DICTIONARY.md` — field and column definitions.
- `FILE_STATUS.csv` — scientific scope and role of each retained file.
- `PROVENANCE.json` — package provenance and proof-status boundary.
- `VALIDATION_REPORT.json` — reproduction and validation summary for this repository update.
- `MANIFEST_SHA256.txt` — SHA-256 manifest for the retained reproducibility package files.

## NPZ layout

`all_responses.npz` contains:

- `N` — `ln(a)`, 1001 samples.
- `k` — 42 wavenumbers in units of `H_d`.
- `values` — shape `(42, 1001, 15)` with columns
  `[zeta_observable, Phi_N, Psi_N, delta_b_N, delta_gamma_N, delta_nu_N, zeta, u, delta_b, delta_gamma, delta_nu, vartheta_b, vartheta_gamma, vartheta_nu, alpha]`.

The first six columns are the observable set used in the stored response maximum. The repeated `zeta` is intentional.

## Scope limits

The package does not provide a CMB likelihood, primordial Boltzmann spectrum, nonlinear galaxy matching, a full matter strong-coupling scale, or a quantum completion. A maximum over the 1001 sampled output times is not a certified continuum maximum. The initial vector is a unit baryon-density response used for linear transfer calculations, not a claim of a physical primordial perturbation of unit amplitude; physical small initial amplitudes rescale the linear response.

The old base-action `S_0` numerical products remain recoverable from `legacy-s0-2026-09-08` and must not be substituted for the curvature-completed observables on `main`.
