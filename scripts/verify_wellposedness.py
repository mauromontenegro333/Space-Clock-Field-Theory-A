#!/usr/bin/env python3
"""Exact algebraic certificate for the SCFT-A nonlinear Cauchy report.

This script checks finite-dimensional identities and exact benchmark
expressions. It deliberately does not claim to machine-prove the analytic
elliptic, pseudodifferential, or quasilinear arguments in the TeX report.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA.mkdir(parents=True, exist_ok=True)


def exact(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(expr)))


def matrix_strings(matrix: sp.Matrix) -> list[list[str]]:
    return [[exact(matrix[i, j]) for j in range(matrix.cols)]
            for i in range(matrix.rows)]


def assert_zero_matrix(name: str, matrix: sp.Matrix) -> None:
    residual = matrix.applyfunc(sp.simplify)
    if residual != sp.zeros(*matrix.shape):
        raise AssertionError(f"{name} failed: {residual}")


def multiset_equal(left: dict[sp.Expr, int], right: dict[sp.Expr, int]) -> bool:
    canonical_left = {sp.factor(sp.simplify(k)): v for k, v in left.items()}
    canonical_right = {sp.factor(sp.simplify(k)): v for k, v in right.items()}
    return canonical_left == canonical_right


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest() -> None:
    ignored_parts = {"build", "preview", "__pycache__"}
    ignored_suffixes = {".zip"}
    records: list[tuple[str, str]] = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(ROOT)
        if any(part in ignored_parts for part in rel.parts):
            continue
        if rel.as_posix() == "data/sha256_manifest.txt":
            continue
        if path.suffix in ignored_suffixes:
            continue
        records.append((sha256(path), rel.as_posix()))
    manifest = "".join(f"{digest}  {name}\n" for digest, name in records)
    (DATA / "sha256_manifest.txt").write_text(manifest, encoding="utf-8")


def main() -> None:
    N, tau, ell, g, rho = sp.symbols(
        "N tau_6 ell g_xi rho", positive=True, finite=True
    )

    T0 = N * sp.Matrix([
        [sp.Rational(8, 3), -sp.Rational(4, 3), 0],
        [-sp.Rational(4, 3), sp.Rational(8, 3), 0],
        [0, 0, 2],
    ])
    V6 = N * sp.Matrix([
        [tau / 2 + ell**2 / (8 * g), tau / 4 + ell**2 / (8 * g), 0],
        [tau / 4 + ell**2 / (8 * g), tau / 2 + ell**2 / (8 * g), 0],
        [0, 0, tau / 2],
    ])
    zero3 = sp.zeros(3)
    A3 = rho**3 * zero3.row_join(T0).col_join((-V6).row_join(zero3))
    S0 = sp.diag(V6, T0)
    principal_residual = sp.simplify(S0 * A3 + A3.T * S0)
    assert_zero_matrix("principal symmetrizer identity", principal_residual)

    frequencies = (T0 * V6).eigenvals()
    expected_frequencies = {
        N**2 * tau: 2,
        N**2 * (tau + ell**2 / (3 * g)): 1,
    }
    if not multiset_equal(frequencies, expected_frequencies):
        raise AssertionError(f"principal frequencies failed: {frequencies}")

    t0_eigenvalues = T0.eigenvals()
    v6_eigenvalues = V6.eigenvals()
    expected_t0 = {sp.Rational(4, 3) * N: 1, 4 * N: 1, 2 * N: 1}
    expected_v6 = {
        N * tau / 4: 1,
        N * (sp.Rational(3, 4) * tau + ell**2 / (4 * g)): 1,
        N * tau / 2: 1,
    }
    if not multiset_equal(t0_eigenvalues, expected_t0):
        raise AssertionError(f"T0 eigenvalues failed: {t0_eigenvalues}")
    if not multiset_equal(v6_eigenvalues, expected_v6):
        raise AssertionError(f"V6 eigenvalues failed: {v6_eigenvalues}")

    # Generic exact Hamiltonian identity with nontrivial grading.
    l11, l12, l22, r1, r2 = sp.symbols(
        "l11 l12 l22 r1 r2", real=True, nonzero=True
    )
    L = sp.Matrix([[l11, l12], [l12, l22]])
    J = sp.Matrix([[0, 1], [-1, 0]])
    R = sp.diag(r1, r2)
    A = R * J * L * R.inv()
    S = R.inv().T * L * R.inv()
    exact_residual = sp.simplify(S * A + A.T * S)
    assert_zero_matrix("exact Hessian symmetrizer identity", exact_residual)

    # Exact inverse of the second-class lapse constraint matrix.
    delta, omega = sp.symbols("Delta_N Omega", nonzero=True)
    constraint_matrix = sp.Matrix([[0, -delta], [delta, omega]])
    claimed_inverse = sp.Matrix([
        [omega / delta**2, 1 / delta],
        [-1 / delta, 0],
    ])
    inverse_left = sp.simplify(constraint_matrix * claimed_inverse)
    inverse_right = sp.simplify(claimed_inverse * constraint_matrix)
    assert_zero_matrix("constraint inverse left residual", inverse_left - sp.eye(2))
    assert_zero_matrix("constraint inverse right residual", inverse_right - sp.eye(2))

    # Repeat with noncommuting exact 2-by-2 blocks to verify operator ordering.
    delta_block = sp.Matrix([[2, 1], [1, 3]])
    omega_block = sp.Matrix([[0, 1], [-1, 0]])
    delta_inverse = delta_block.inv()
    zero2 = sp.zeros(2)
    block_constraint = zero2.row_join(-delta_block).col_join(
        delta_block.row_join(omega_block)
    )
    block_inverse = (
        (delta_inverse * omega_block * delta_inverse).row_join(delta_inverse)
        .col_join((-delta_inverse).row_join(zero2))
    )
    noncommuting_left = sp.simplify(block_constraint * block_inverse)
    noncommuting_right = sp.simplify(block_inverse * block_constraint)
    assert_zero_matrix(
        "noncommuting constraint inverse left residual",
        noncommuting_left - sp.eye(4),
    )
    assert_zero_matrix(
        "noncommuting constraint inverse right residual",
        noncommuting_right - sp.eye(4),
    )

    # Principal lapse Schur sign.
    r_scalar = sp.symbols("r_xi", real=True)
    B4 = -N * ell * rho**2 * r_scalar / 2
    C2 = -2 * N * g * rho**2
    schur_term = sp.factor(-B4**2 / C2)
    expected_schur = N * ell**2 * rho**2 * r_scalar**2 / (8 * g)
    if sp.simplify(schur_term - expected_schur) != 0:
        raise AssertionError("lapse Schur-complement sign failed")

    # Exact manuscript benchmark; all expressions remain symbolic/rational.
    Q_star = sp.Rational(101, 100)
    N_star = 1 / Q_star
    U2_over_Hd2 = sp.Integer(10)**36 / (
        sp.Integer(9999) * sp.sqrt(sp.Integer(9999))
    )
    m2_over_Hd2 = sp.factor(Q_star**2 * U2_over_Hd2)
    expected_m2 = (
        sp.Integer(1020100) * sp.Integer(10)**30
        / (sp.Integer(9999) * sp.sqrt(sp.Integer(9999)))
    )
    if sp.simplify(m2_over_Hd2 - expected_m2) != 0:
        raise AssertionError("exact rolling mass identity failed")
    if not (N_star > 0 and U2_over_Hd2 > 0 and m2_over_Hd2 > 0):
        raise AssertionError("benchmark positivity failed")

    symbol_payload = {
        "assumptions": [
            "N > 0", "tau_6 > 0", "ell > 0", "g_xi > 0", "rho > 0"
        ],
        "T0": matrix_strings(T0),
        "V6": matrix_strings(V6),
        "T0_eigenvalues": {
            exact(key): value for key, value in t0_eigenvalues.items()
        },
        "V6_eigenvalues": {
            exact(key): value for key, value in v6_eigenvalues.items()
        },
        "squared_frequencies": {
            exact(key): value for key, value in frequencies.items()
        },
        "principal_symmetrizer_residual": matrix_strings(principal_residual),
        "exact_hessian_symmetrizer_residual": matrix_strings(exact_residual),
        "constraint_matrix_left_residual": matrix_strings(
            inverse_left - sp.eye(2)
        ),
        "constraint_matrix_right_residual": matrix_strings(
            inverse_right - sp.eye(2)
        ),
        "noncommuting_constraint_left_residual": matrix_strings(
            noncommuting_left - sp.eye(4)
        ),
        "noncommuting_constraint_right_residual": matrix_strings(
            noncommuting_right - sp.eye(4)
        ),
        "positive_lapse_schur_term": exact(schur_term),
    }
    write_json(DATA / "exact_symbol_checks.json", symbol_payload)

    benchmark_rows = [
        ("Q_star", exact(Q_star), "dimensionless", "positive"),
        ("N_star", exact(N_star), "dimensionless", "positive"),
        ("nu_star_Hd_squared", "1/(2*10^170)", "dimensionless", "positive"),
        ("ell_star_Hd_squared", "1/(2*10^30)", "dimensionless", "positive"),
        ("tau_6_Hd_fourth", "10^(-340)", "dimensionless", "positive"),
        ("U_star_second_over_Hd_squared", exact(U2_over_Hd2),
         "dimensionless", "positive"),
        ("m_star_squared_over_Hd_squared", exact(m2_over_Hd2),
         "dimensionless", "positive"),
    ]
    with (DATA / "exact_benchmark.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "exact_expression", "units", "sign"])
        writer.writerows(benchmark_rows)

    certificate = {
        "certificate_version": 1,
        "date_utc": "2026-09-09",
        "theorem_scope": {
            "sector": "vacuum gravity-clock sector of the completed SCFT-A action",
            "spatial_domain": "compact boundaryless three-torus",
            "data_domain": "open strongly Dirac-regular component",
            "topology": "H^(s+3) x H^s for integer s >= 18",
            "result": [
                "local existence",
                "uniqueness",
                "continuous dependence",
                "persistence of higher regularity",
                "constraint propagation",
                "continuation criterion",
            ],
            "excluded_claims": [
                "global existence",
                "long-time nonlinear stability",
                "bounded-domain initial-boundary-value well-posedness",
                "unconditional evolution for unspecified S_vis",
                "quantum or strong-coupling control",
            ],
        },
        "exact_checks": {
            "principal_symmetrizer_identity": "PASS",
            "principal_frequency_multiplicities": "PASS",
            "momentum_hessian_positive_eigenvalues": "PASS",
            "metric_hessian_positive_eigenvalues": "PASS",
            "exact_hessian_symmetrizer_identity": "PASS",
            "second_class_constraint_inverse": "PASS",
            "lapse_schur_sign": "PASS",
            "rolling_benchmark_identity_and_signs": "PASS",
        },
        "machine_check_boundary": (
            "The script certifies exact algebra and artifact integrity only. "
            "The functional-analytic proof is the TeX/PDF report and remains "
            "subject to expert mathematical review."
        ),
    }
    write_json(DATA / "proof_certificate.json", certificate)
    write_manifest()
    print("PASS: 8 exact algebraic/structural checks")
    print(f"Wrote {DATA / 'exact_symbol_checks.json'}")
    print(f"Wrote {DATA / 'exact_benchmark.csv'}")
    print(f"Wrote {DATA / 'proof_certificate.json'}")
    print(f"Wrote {DATA / 'sha256_manifest.txt'}")


if __name__ == "__main__":
    main()
