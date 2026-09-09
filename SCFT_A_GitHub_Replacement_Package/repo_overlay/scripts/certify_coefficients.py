#!/usr/bin/env python3
"""Rigorous sufficient bounds, not an approximation to the old extremum.

All inequalities except the tilted-mean enclosure use fractions.Fraction.
The one integral enclosure uses Arb outward-rounded range Riemann sums.
No floating-point optimizer, quadrature error estimate, or sampled minimum
is used in the certificate. See the revised technical companion for the
analytic bounds between and outside these integration boxes.
"""
from fractions import Fraction as F
from pathlib import Path
import json
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]

def tilted_mean_enclosure(n=4096, bits=160):
    ctx.prec = bits
    I0, I1 = arb(0), arb(0)
    for i in range(n):
        lo, hi = arb(i)/n, arb(i+1)/n
        z = lo.union(hi)
        if i == 0:
            # b increases on this box; exp(-16 z) <= 1.
            p = arb(0).union((-1/(hi*(1-hi))).exp())
        elif i == n-1:
            # Both b and exp(-16 z) decrease on this box.
            p = arb(0).union((-1/(lo*(1-lo))-16*lo).exp())
        else:
            p = (-1/z-1/(1-z)-16*z).exp()
        I0 += p/n
        I1 += z*p/n
    mean = I1/I0
    assert I0 > 0 and mean > arb(1)/4 and mean < arb(3)/10
    return {"boxes": n, "precision_bits": bits,
            "I0_ball": str(I0), "I1_ball": str(I1),
            "mean_ball": str(mean), "certified_mean_interval": ["1/4", "3/10"]}

def main():
    Qa, Qb, L, eps = F(1001,1000), F(101,100), F(9,1000), F(1,10**6)
    checks = {}
    def check(name, condition):
        checks[name] = bool(condition)
        assert condition, name
    # Moment target, using exact rational enclosures on Fb and ub.
    zbar_lo = 1 - F(1,19900)/(L*(F(1,100)-2*eps))
    zbar_hi = 1 - (F(1,20000)-2*eps**2)/(L*F(1,99))
    check("target_mean_lower", zbar_lo > F(2,5))
    check("target_mean_upper", zbar_hi < F(1,2))
    check("bump_normalization_lower", (F(1,100)-2*eps)/L > 1)
    check("bump_normalization_upper", 2*3**18/(99*L) < 10**9)
    # On z <= 1/16, R=(log phi)' >= (9/10) z^-2;
    # -R' <= (201/100) z^-3, hence -R'/R^2 <= (5/2)z.
    check("inner_log_slope", 1-F(18,256) > F(9,10))
    check("inner_log_curvature", 2+2/F(15,1)**3 < F(201,100))
    check("inner_ratio", F(201,100)/F(9,10)**2 < F(5,2))
    check("inner_am_margin", F(96,37) > 1+F(1515,1001))
    hmin = F(1,10**22)
    check("compact_h_lower", Qa*L/(96*3**35) > hmin)
    common = Qa*L/(3*F(37,32))
    check("inner_d0_lower", 8*common*256/(Qb**4*10**9) > F(1,10**12))
    check("inner_d1_lower", common/(Qb**3*10**5) > F(1,10**12))
    check("compact_bump_log_slope", (16**2+18000**2+16)/L < 4*10**10)
    check("compact_ramp_log_slope", 1+(20000**2+4)/eps < 5*10**14)
    check("ramp_start_exponential_bound", 3*10**6*(20000**2+4) < 2**10978)
    check("global_g_derivative", 18*10**9/L + 10**7 < 3*10**12)
    hmax, mmax, mnmax, rmax = F(1,250), F(2*10**9), F(10**14), F(10**10)
    qmax = F(1515,1001)
    check("global_h_upper", Qb/F(297) < hmax)
    check("global_m_upper", Qb**2*(F(10**9,16)+2) < mmax)
    check("global_mN_upper", 6*Qb/99+3*Qb**2*10**15/99 < mnmax)
    check("rN_upper", F(20,1)/F(3,1000)*F(1,99)*3**16 < rmax)
    A1 = mmax*(qmax-1)+mnmax+rmax*mmax
    E2 = (6*hmax*(mmax*(qmax+1)+mnmax)+mmax**2
          +4*hmax*(mmax*(qmax-1)+mnmax)+4*hmax*rmax*mmax+8*hmax*rmax)
    E3 = 6*hmax*(mmax*(qmax+1)+mnmax)+2*mmax+8*hmax*rmax
    check("A1_upper", A1 < 10**20)
    check("E2_upper", E2 < 10**20)
    check("E3_upper", E3 < 10**14)
    delta = F(1,10**140)
    cmin, dmin = 2*hmin**2, hmin
    # Stronger than the original 1/2 min prescription, without finding extrema.
    check("delta_unit", 2*delta <= 1)
    check("delta_inner_d0", 2*delta <= F(1,10**12))
    check("delta_inner_d1", 2*delta <= F(1,10**12))
    check("delta_A1", 2*delta*10**20 <= hmin/2)
    check("delta_E2", (2*delta)**2 <= 27*cmin**2*dmin/(32*F(10**20)**3))
    check("delta_E3", 2*delta <= 27*cmin*dmin**2/(32*F(10**14)**3))
    # Explicit regulator: uniform tip bound <= 40000 sigma.
    sigma = F(1,10**6)
    check("tip_size", 77+10*3**6*5 < 40000)
    check("tip_hessian", 1-40000*sigma > F(1,2))
    Ck = F(10**30)
    check("endpoint_kinetic", 500*Ck-12*delta > 0)
    check("endpoint_energy", Ck*delta < F(1,10))
    mean = tilted_mean_enclosure()
    checks["tilted_mean_interval"] = True
    report = {"all_passed": all(checks.values()), "checks": checks,
              "tilted_mean_certificate": mean,
              "exact_parameters": {"kappa_m_over_Hd2": "10^30", "delta": "1/10^140",
                 "sigma": "1/10^6", "s_N": "10", "a0_over_Hd": "1"},
              "certified_bounds": {"theta": "(-16,0)", "c": "(1,10^9)",
                 "n0": "16", "h_star_lower": "1/10^22", "d0_lower": "1/10^12",
                 "d1_lower": "1/10^12", "A1_upper": "10^20", "E2_upper": "10^20", "E3_upper": "10^14"},
              "scope": "Sufficient bounds for a NEW exact coefficient point; not certification of the old numerical extremum."}
    (ROOT/"data").mkdir(exist_ok=True)
    (ROOT/"data/coefficient_certificate.json").write_text(json.dumps(report,indent=2)+"\n")
    print(f"{len(checks)} coefficient checks passed")

if __name__ == "__main__":
    main()
