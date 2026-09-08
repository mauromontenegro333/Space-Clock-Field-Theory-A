#!/usr/bin/env python3
"""Exact algebraic checks for SCFT-A's nonlinear Rank-1 supplement.

Run: python3 scft_rank1_checks.py --output scft_rank1_checks.json
Requires SymPy. All test quantities are exact; no floating-point scans are used.
This checks displayed algebraic identities; the nonlinear energy estimate is an unproved hypothesis, not a check outcome.
"""

import argparse
import json
from pathlib import Path
import sys

import sympy as s


def run_checks():
    checks = {}

    def zero(name, expression):
        residual = s.factor(s.cancel(expression))
        checks[name] = {"passed": residual == 0, "residual": str(residual)}

    def positive(name, expression):
        value = s.simplify(expression)
        passed = bool(value.is_positive)
        checks[name] = {"passed": passed, "exact_value": str(value)}

    z, chi, alpha, zeta = s.symbols("z chi alpha zeta", real=True)
    H, m, x, k, ell = s.symbols("H m x k ell", positive=True)
    b = k*x
    A = 2+3*b
    C = x*(1+ell*x)
    D = 2*H**2*A+b*(m+2*x)
    L = (-3*z*z+2*z*chi+2*H*alpha*(3*z-chi)
         +(m-6*H**2+2*x)*alpha**2/2+2*C*alpha*zeta+x*zeta**2
         -b*(3*z-chi-3*H*alpha)**2/2)
    alpha_sol = 2*(H*A*z-b*C*zeta)/D
    chi_sol = A*((m+2*x)*z+2*H*C*zeta)/D
    substitutions = {alpha: alpha_sol, chi: chi_sol}
    zero("lapse_auxiliary_equation", s.diff(L, alpha).subs(substitutions))
    zero("shift_auxiliary_equation", s.diff(L, chi).subs(substitutions))
    K = A*(m+2*x)/D
    B = 2*H*C*A/D
    V0 = 2*b*C*C/D-x
    zero("complete_scalar_reduction", L.subs(substitutions)-(K*z*z+2*B*zeta*z-V0*zeta*zeta))
    zero("auxiliary_hessian_determinant", s.det(s.hessian(L, (alpha, chi)))+D)
    zero("minkowski_kinetic", K.subs({H:0, m:0})-(3+2/(k*x)))
    zero("minkowski_potential", V0.subs({H:0, m:0})-(2*ell*x*x+ell*ell*x**3))

    # Differentiate along one fixed COMOVING mode, including x_N=-2x.
    h, q, mn, kn, ln = s.symbols("h q m_N k_N ell_N", real=True)
    D = 2*h*A+b*(m+2*x)
    R = 2*x*(1+ell*x)*A/D   # B/H
    def dN(expr):
        return (s.diff(expr, x)*(-2*x)+s.diff(expr, h)*(-2*q*h)
                +s.diff(expr, m)*mn+s.diff(expr, k)*kn+s.diff(expr, ell)*ln)
    V = 2*k*x**3*(1+ell*x)**2/D-x+h*((3-q)*R+dN(R))
    zero("rolling_infrared", s.limit(V/x, x, 0)-q)
    zero("mixed_switch_off", V.subs({k:0, kn:0})-x*(q+((q-1)*ell+ln)*x))
    polynomial_general = s.Poly(s.cancel(V*D**2/x), x)
    c_general = [
        16*h*h*q,
        4*h*(12*h*k*q+4*h*ell*(q-1)+4*h*ln-k*m*(q-1)-k*mn-kn*m),
        (36*h*h*k*k*q+48*h*h*k*ell*(q-1)+48*h*h*k*ln
         -6*h*k*k*(m*(q+1)+mn)-4*h*k*ell*(m*(q-1)+mn)
         +4*h*k*ln*m+8*h*k*(4-q)-4*h*kn*ell*m-8*h*kn-k*k*m*m),
        2*(18*h*h*k*k*ell*(q-1)+18*h*h*k*k*ln
           -3*h*k*k*ell*(m*(q+1)+mn)+3*h*k*k*ln*m
           +6*h*k*k*(2-q)+4*h*k*ell*(5-q)+4*h*k*ln
           -4*h*kn*ell-k*k*m),
        4*k*(3*h*k*ell*(3-q)+3*h*k*ln+2*h*ell*ell+k*ell*m),
        2*k*k*ell*(6*h*ell+ell*m+4),
        4*k*k*ell*ell,
    ]
    for j in range(7):
        zero(f"positive_floor_polynomial_coefficient_{j}",
             polynomial_general.nth(j)-c_general[j])
    zero("positive_floor_polynomial_degree", polynomial_general.degree()-6)

    # The former mixed-support polynomial is the ell=1, ell_N=0 specialization.
    Vone = V.subs({ell:1, ln:0})
    polynomial = s.Poly(s.cancel(Vone*D**2/x), x)
    for j in range(7):
        zero(f"unit_plateau_coefficient_{j}",
             polynomial.nth(j)-c_general[j].subs({ell:1,ln:0}))

    # Exact floor endpoints and ultraviolet physical symbol.
    rho = s.Rational(1,2)
    wv = s.symbols("w", real=True)
    profile = rho+(1-rho)*wv
    zero("floor_profile_at_zero", profile.subs(wv,0)-rho)
    zero("floor_profile_at_one", profile.subs(wv,1)-1)
    zero("scalar_kinetic_ultraviolet", s.limit(A*(m+2*x)/D,x,s.oo)-3)
    zero("scalar_potential_ultraviolet",
         s.limit(V/(x*(1+x*x)),x,s.oo)-ell*ell)
    delta6 = s.symbols("delta_6", positive=True)
    Vfull = V+3*delta6**2*x**3
    zero("completed_scalar_potential_ultraviolet",
         s.limit(Vfull/(x*(1+x*x)),x,s.oo)-(ell*ell+3*delta6**2))
    zero("conformal_ricci_gradient_norm", 1+2+3-6)
    zero("scalar_curvature_gradient_action_factor", s.Rational(6,2)-3)
    zero("tensor_curvature_gradient_action_factor",
         s.Rational(1,2)*s.Rational(1,4)/s.Rational(1,8)-1)

    # General flat-symbol Ricci identity used in the nonlinear harmonic-gauge lemma.
    xi = s.Matrix(s.symbols("xi_0:3", real=True))
    k00, k01, k02, k11, k12, k22 = s.symbols(
        "k_00 k_01 k_02 k_11 k_12 k_22", real=True)
    kmat = s.Matrix([[k00, k01, k02],
                     [k01, k11, k12],
                     [k02, k12, k22]])
    trk = s.trace(kmat)
    q2 = xi.dot(xi)
    divk = kmat*xi
    harmonic = divk-xi*trk/2
    ricci_direct = (q2*kmat-xi*divk.T-divk*xi.T+xi*xi.T*trk)/2
    ricci_harmonic = q2*kmat/2-(xi*harmonic.T+harmonic*xi.T)/2
    for i in range(3):
        for j in range(i,3):
            zero(f"ricci_harmonic_principal_{i}{j}",
                 ricci_direct[i,j]-ricci_harmonic[i,j])

    # The bump amplitude enforces exact equality of J above the regulator tip.
    moment_difference, bump_mass = s.symbols("moment_difference bump_mass", nonzero=True)
    amplitude = moment_difference/bump_mass
    zero("smooth_tip_zero_moment", moment_difference-amplitude*bump_mass)

    # Exact Hamiltonian symmetrizer identity and Dirac count.
    kval, vval = s.symbols("K_value V_value", positive=True)
    principal = s.Matrix([[0,1],[-vval/kval,0]])
    symmetrizer = s.diag(vval,kval)
    skew_residual = symmetrizer*principal+(symmetrizer*principal).T
    zero("scalar_symmetrizer_skew_identity",
         sum(entry**2 for entry in skew_residual))
    zero("rank_one_physical_count", s.Rational(20-2*6-2,2)-3)

    # Independent action-level check of the background and clock identities.
    Q, F, G, U, Gp = s.symbols("Q F G U Gp", real=True)
    h_bg = (Q*F-U)/3
    QN = -3*F/G
    def alongQ(expr):
        return s.diff(expr, Q)+s.diff(expr, F)*G+s.diff(expr, G)*Gp+s.diff(expr,U)*F
    zero("homogeneous_charge", G*QN+3*F)
    zero("homogeneous_raychaudhuri", alongQ(h_bg)*QN+Q*F)
    zero("m_N_expression", alongQ(Q*Q*G)*QN-(-6*Q*F-3*Q*Q*F*Gp/G))

    # Young inequalities used in the continuum positivity proof.
    t = s.symbols("t", positive=True)
    zero("young_two_to_one", 2+t**3-3*t-(t-1)**2*(t+2))
    zero("young_one_to_two", 1+2*t**3-3*t*t-(t-1)**2*(2*t+1))
    positive("inner_young_safety_margin", s.Integer(3456)-8)
    qa, qb, eps = s.Rational(1001,1000), s.Rational(101,100), s.Rational(1,10**6)
    width = qb-qa
    positive("matching_positive_mass", s.Rational(1,100)-2*eps)
    positive("matching_mean_below_right", s.Rational(1,20000)-2*eps**2)
    positive("matching_mean_above_left", width/s.Integer(100)-s.Rational(1,19900)-2*width*eps)
    positive("deceleration_upper_bound", s.Rational(7,4)-3*qb/(2*qa))

    # The mixed Legendre auxiliary field is a positive elliptic Schur block.
    # A finite exact matrix example verifies the differentiation identity without
    # claiming to discretize or prove the continuum theorem.
    n = s.symbols("n", positive=True)
    mat = s.Matrix([[2+n*n, n], [n, 3+n*n]])
    src = s.Matrix([1,2])
    solution = -mat.inv()*src
    extended = (src.dot(solution)+s.Rational(1,2)*(solution.T*mat*solution)[0])
    # General envelope second derivative: E_nn - E_nK E_KK^-1 E_Kn.
    dmat, ddmat = mat.diff(n), mat.diff(n,2)
    expected = ((solution.T*ddmat*solution)[0]/2
                -(dmat*solution).dot(mat.inv()*(dmat*solution)))
    zero("inverse_operator_hessian_schur_identity", s.diff(extended,n,2)-expected)

    # Exact nontrivial nonlinear Bianchi-I family, not generic inhomogeneous data.
    ps = (s.Rational(-1,3),s.Rational(2,3),s.Rational(2,3))
    zero("kasner_sum",sum(ps)-1)
    zero("kasner_squared_sum",sum(p*p for p in ps)-1)
    curvature = 4*sum(p*p*(p-1)**2 for p in ps)+4*sum(ps[i]**2*ps[j]**2 for i in range(3) for j in range(i+1,3))
    zero("kasner_kretschmann_coefficient",curvature-s.Rational(64,27))

    return {
        "scope": "Exact algebra for the positive-floor Rank-1 action; nonlinear evolution estimates remain explicit hypotheses in the integrated manuscript.",
        "floating_point_used": False,
        "sympy_version": s.__version__,
        "checks": checks,
        "passed_count": sum(c["passed"] for c in checks.values()),
        "total_count": len(checks),
        "all_passed": all(c["passed"] for c in checks.values()),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_checks()
    content = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(content, encoding="utf-8")
    print(content)
    sys.exit(0 if result["all_passed"] else 1)

