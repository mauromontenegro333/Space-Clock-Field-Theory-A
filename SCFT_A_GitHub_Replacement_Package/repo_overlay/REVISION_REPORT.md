# Current status of the eight major SCFT-A issues

This report combines the revised classical-consistency package with the later local nonlinear vacuum Cauchy result. The manuscript files themselves are intentionally not distributed in this repository.

| # | Issue | Current status | Remaining boundary |
|---:|---|---|---|
| 1 | Full nonlinear well-posedness | **Locally established in the stated vacuum domain.** The newer theorem gives local existence, uniqueness, continuous dependence, persistence of higher regularity, constraint propagation, and a continuation criterion on an open strongly Dirac-regular component of a compact boundaryless three-torus, in `H^(s+3) x H^s` for integer `s >= 18`. | No global existence/long-time nonlinear stability theorem; no bounded-domain initial-boundary-value theorem; no unconditional extension to unspecified `S_vis`. |
| 2 | EFT and interaction control | The specified classical ideal-fluid cubic density and gravitational cubic functional are available in the revised work. | No complete strong-coupling domain, complete quartic/exchange analysis, radiative protection, or quantum completion. |
| 3 | Stability beyond instantaneous signs | An exact decreasing endpoint scalar energy is established on the stated expanding de Sitter endpoint domain. | General time-dependent-background and global nonlinear stability are not established. |
| 4 | AQUAL/MOND coefficient hierarchy | Repaired at the exact scale-separated point `kappa_m/H_d^2 = 10^30`, with the declared stationary correction bound and recomputed backgrounds/responses. | Nonlinear sourced cosmological-to-galactic matching remains open. |
| 5 | Contracting reciprocal branch | Explicit source-free unbraided contracting branch supplied, with `H < 0` and the stated clock-rate behavior. | No contracting extension through nonzero braiding is claimed. |
| 6 | Claim discipline | Classical results, theorem domains, numerical approximations, and open problems are separated explicitly. | No completed dark-sector-unification or publication guarantee is claimed. |
| 7 | Engineered functions and coefficients | Plateau, switches, positive floors, curvature term, and scale separation are identified as a deliberate mathematical existence construction. | No microscopic/naturalness mechanism selects the hierarchy. |
| 8 | Uncertified extrema | The old inherited extremum is not used. The revised point uses exact `delta = 10^-140` and `sigma = 10^-6` with sufficient-bound certificates. | The discarded old extremum is not retroactively certified. |

The nonlinear theorem certificate is under `nonlinear_wellposedness/data/proof_certificate.json`. Exact finite-dimensional symbol identities are under `nonlinear_wellposedness/data/exact_symbol_checks.json`. The theorem remains an analytic mathematical result subject to expert review; the verifier is not a proof assistant for the functional-analysis argument.
