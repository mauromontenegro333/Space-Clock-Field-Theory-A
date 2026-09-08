# SCFT integrated manuscript data dictionary

The numerical CSV and NPZ files in this package support the retained base-action $S_0$ calculations. They do not describe finite-momentum perturbations of the completed spatial action. The completed-action material in `checks/` consists of exact algebraic identity reports.

All array columns below are **zero-indexed**. Each stored value is a binary64 numerical approximation. Raw delivered file contents are preserved. The numerical normalization uses `H_d=1`, densities divided by `M_Pl^2 H_d^2`, and conserved charge normalization `C_T/H_d^2=1`. It does not set the physical ratio `M_Pl/H_d` to one. `N = ln(a)`, and primes are derivatives with respect to `N`.

## Background CSV files

The loaded background files have a comment header beginning with `#`. The source-free CSV has no header. All have the same 24-column layout and are readable with `np.loadtxt(path, delimiter=',')`.

| Column | Header | Meaning |
| --- | --- | --- |
| 0 | `N` | `ln(a)`. |
| 1 | `Q` | Dimensionless homogeneous clock norm. |
| 2 | `H` | `H/H_d`. |
| 3 | `Qprime` | `dQ/dN`. |
| 4 | `Hprime` | `d(H/H_d)/dN`. |
| 5 | `Xi` | `Xi_s/H_d^2`, where `Xi_s = D + 3 B_Q^2/2`. |
| 6 | `W` | `W_s/H_d`, where `W_s = 2H - Q B_Q`. |
| 7 | `C` | `C_s/H_d^2`, the scalar gradient numerator. |
| 8 | `Sigma` | Acceleration regulator `Sigma_R(Q)`. |
| 9 | `Upsilon` | `1 - 2 Sigma_R (rho_b + 4 rho_r/3)/(M_Pl^2 W_s^2)`. |
| 10 | `kmix_over_aH` | Positive density-reduction boundary `k_mix/(aH)`, or zero when the implemented radicand is nonpositive. |
| 11 | `D0_over_H2` | `Q^2 Xi_s/H^2`. |
| 12 | `W_over_H` | `W_s/H`. |
| 13 | `A_over_H2` | `[Q^2 D - 6H^2 + 6HQ B_Q]/H^2`; here `D = K_QQ + 3H B_QQ`. |
| 14 | `C_over_H2` | `C_s/H^2`. |
| 15 | `rhob_over_H2` | `rho_b/(M_Pl^2 H^2)`. |
| 16 | `rhor_over_H2` | `rho_r/(M_Pl^2 H^2)`. |
| 17 | `K` | `K_R/H_d^2`. |
| 18 | `K1` | `K_R,Q/H_d^2`; this equals the conserved charge density only on the unbraided branch. |
| 19 | `K2` | `K_R,QQ/H_d^2`. |
| 20 | `K3` | `K_R,QQQ/H_d^2`. |
| 21 | `K4` | `K_R,QQQQ/H_d^2`. |
| 22 | `Sigma1` | `d Sigma_R/dQ`. |
| 23 | `Sigma2` | `d^2 Sigma_R/dQ^2`. |

The source-free background has visible densities zero and `Upsilon=1`. Its time interval is longer than the loaded runs. `data/reference/scft_background.csv` has 20,001 rows and is shared by the reference and scan results; the exact duplicate scan copy was omitted. `data/refinement/scft_background.csv` has 40,001 rows; the source-free CSV has 4,001 rows.

## Mode NPZ files

Each NPZ is a compressed dictionary of NumPy arrays. Keys are the solver's decimal strings for `k/H_d`. Load with `allow_pickle=False`. Each mode array has 6,001 rows on a uniform `N` grid from `ln(1e-10)` through `ln(10)`.

| Column | Quantity | Meaning |
| --- | --- | --- |
| 0 | `N` | `ln(a)`. |
| 1 | `zeta` | Unitary-gauge clock curvature divided by the initial normalization `zeta_*`. |
| 2 | `Psi_N` | Newtonian-gauge curvature potential divided by `zeta_*`. |
| 3 | `delta_b_N` | Newtonian-gauge baryon density contrast divided by `zeta_*`. |
| 4 | `delta_r_N` | Newtonian-gauge radiation density contrast divided by `zeta_*`. |
| 5 | `P_over_H2` | Density-only reduction factor `mathscr P_k/H^2`; it is the unbounded quantity before plotting `P/(1+abs(P))`. |
| 6 | `k_over_aH` | Physical wavenumber in Hubble units, `k/(aH)`. |

The integrations set `zeta_*=1` to compute transfer functions. This does not assign unit-amplitude physical cosmological perturbations. The mode archives store the seven reported quantities; they do not store every internal canonical variable or the integrator's adaptive steps.

## Main audit JSON

| Key | Meaning |
| --- | --- |
| `FE`, `KE` | Numerical matching values `K_Q(Q_E)/H_d^2` and `K(Q_E)/H_d^2`. |
| `background_samples` | Number of background grid rows used for interpolation. |
| `minimum_Upsilon`, `maximum_kmix_aH`, `minimum_Xi`, `minimum_W`, `minimum_C` | Extrema on that stored background sample grid. |
| `maximum_Qprime`, `maximum_Hprime` | Largest sampled derivatives, useful for checking the signs on the sampled interval. |
| `modes` | List of per-mode summaries described below. |

| Per-mode key | Meaning |
| --- | --- |
| `k` | `k/H_d`. |
| `max_r` | Maximum `k/(aH)` on the 6,001-point output grid. |
| `crossings` | Number of adjacent output samples with opposite signs of `mathscr P_k`; not a continuum root certificate. |
| `db1`, `db10` | Baryon transfer at `a=1` and `a=10`. |
| `psi1` | Newtonian curvature transfer at `a=1`. |
| `maxamp` | Largest absolute value among `zeta`, `Psi_N`, `delta_b_N`, and `delta_r_N` on the output grid. |
| `nfev` | Integrator-reported right-hand-side evaluations. |

`refinement_differences.json` uses, for each listed quantity, `max_i(abs(fine_i-base_i)/max(1,abs(fine_i)))`. Mode `10.0` uses the reference baseline; mode `362.0291470135738` uses the scan baseline.

`regulator_checks.json` records the optimized values of `abs(Q Sigma_Q/Sigma)`, `abs(Q^2 Sigma_QQ/Sigma)`, their locations, `max(abs(Sigma_dot/(H Sigma)))`, and `min(1/G1,1/sqrt(G2))`. `refined_background.json` contains refined rather than grid-only extrema; its values can therefore differ slightly from those in the main audit JSON. Equality and acceleration locations use the manuscript's internal scale-factor normalization.

`symbolic_checks.json` is the retained Boolean result record. The corresponding interactive command text was not preserved as a standalone script; it must not be mistaken for an executable proof.

## Completed-action exact check reports

`checks/scft_completion_checks.json` records 14 exact completed auxiliary and quadratic identities. `checks/scft_rank1_checks.json` records 54 exact local-transition identities and positivity calculations. Both reports state `floating_point_used: false`. These algebraic checks do not establish the loss-free nonlinear energy estimate or global constraint compatibility required for full nonlinear Rank-1 closure.
