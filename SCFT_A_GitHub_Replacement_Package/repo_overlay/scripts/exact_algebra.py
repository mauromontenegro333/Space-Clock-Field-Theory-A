#!/usr/bin/env python3
"""Independent symbolic elimination and positivity/energy identities."""
from pathlib import Path
import json
import sympy as s
ROOT = Path(__file__).resolve().parents[1]

def main():
    checks={}
    def zero(name, expr):
        ok=s.cancel(expr)==0
        checks[name]=bool(ok)
        assert ok,name
    z,al,chi,Z,H,W,ak,b,C,Rv=s.symbols('z al chi Z H W ak b C Rv')
    A=2+3*b;J=W+3*b*H;D=W**2+b*(ak+6*H*W)
    L=-3*z**2+W*al*(3*z-chi)+2*z*chi+ak*al**2/2+2*C*al*Z+Rv*Z**2-b*(3*z-chi-3*H*al)**2/2
    K=(A*ak+3*W**2+18*b*H*(W-H))/D
    sol={al:(2*J*z-2*b*C*Z)/D,chi:K*z+2*J*C*Z/D}
    zero('lapse_constraint',s.diff(L,al).subs(sol))
    zero('shift_constraint',s.diff(L,chi).subs(sol))
    zero('auxiliary_determinant',s.det(s.hessian(L,[al,chi]))+D)
    zero('reduced_action',L.subs(sol)-(K*z**2+4*J*C/D*z*Z+(Rv-2*b*C**2/D)*Z**2))
    x,d,c=s.symbols('x delta c',positive=True)
    bd=d*x/(2*c);Ad=2+3*bd;Wd=s.Rational(6,5);akd=2*x-s.Rational(158,125)
    Dd=Wd**2+bd*(akd+6*Wd);Cd=x*(1+x/(2*c));Jd=Wd+3*bd
    Kd=s.factor((Ad*akd+3*Wd**2+18*bd*(Wd-1))/Dd)
    Bd=2*Jd*Cd/Dd
    Vd=2*bd*Cd**2/Dd-x+4*x**2/c+3*d**2*x**3/c**2+3*Bd-2*x*s.diff(Bd,x)
    Zn=180*c+125*d*x**2+371*d*x
    zero('endpoint_denominator',Dd-Zn/(125*c))
    zero('endpoint_kinetic',Kd-(500*c*x+224*c+375*d*x**2-12*d*x)/Zn)
    nV,denV=s.fraction(s.factor(Vd))
    checks['endpoint_positive_polynomial']=all(a>0 for a in s.Poly(nV,x,d,c).coeffs())
    assert checks['endpoint_positive_polynomial']
    zero('endpoint_small_delta',s.limit(Vd,d,0)-(s.Rational(2,3)*x+s.Rational(19,6)*x**2/c))
    nk,dk=s.fraction(s.factor(3*Kd-x*s.diff(Kd,x)))
    nv,dv=s.fraction(s.factor(s.diff(Vd,x)))
    # Prove energy inequalities at the exact revised point by exact coefficient signs.
    pars={c:10**30,d:s.Rational(1,10**140)}
    for name,poly in [('endpoint_kinetic_dissipation',nk),('endpoint_restoring_monotonicity',nv)]:
        checks[name]=all(co>0 for co in s.Poly(poly.subs(pars),x).all_coeffs())
        assert checks[name],name
    # Independently recover the local-transition polynomial with both switch derivatives.
    h,m,q,k,lam,mn,kn,ln=s.symbols('h m q k lam mn kn ln',positive=True)
    aa=2+3*k*x;dd=2*h*aa+k*x*(m+2*x);cc=x*(1+lam*x)
    bb=2*s.sqrt(h)*cc*aa/dd
    DN=lambda f:-2*q*h*s.diff(f,h)-2*x*s.diff(f,x)+mn*s.diff(f,m)+kn*s.diff(f,k)+ln*s.diff(f,lam)
    vv=2*k*x*cc**2/dd-x+s.sqrt(h)*(DN(bb)+3*bb)
    co=[16*h**2*q,
      4*h*(12*h*k*q+4*h*lam*(q-1)+4*h*ln-k*m*(q-1)-k*mn-kn*m),
      36*h**2*k**2*q+48*h**2*k*lam*(q-1)+48*h**2*k*ln-6*h*k**2*(m*(q+1)+mn)-4*h*k*lam*(m*(q-1)+mn)+4*h*k*ln*m+8*h*k*(4-q)-4*h*kn*lam*m-8*h*kn-k**2*m**2,
      2*(18*h**2*k**2*lam*(q-1)+18*h**2*k**2*ln-3*h*k**2*lam*(m*(q+1)+mn)+3*h*k**2*ln*m+6*h*k**2*(2-q)+4*h*k*lam*(5-q)+4*h*k*ln-4*h*kn*lam-k**2*m),
      4*k*(3*h*k*lam*(3-q)+3*h*k*ln+2*h*lam**2+k*lam*m),
      2*k**2*lam*(6*h*lam+lam*m+4),4*k**2*lam**2]
    zero('local_polynomial_both_switches',dd**2*vv/x-sum(a*x**i for i,a in enumerate(co)))
    Q,F,g,U=s.symbols('Q F g U',positive=True)
    Qdot=-3*H*F/g;Hdot=-Q*F/2
    zero('unbraided_charge_propagation',g*Qdot+3*H*F)
    zero('unbraided_friedmann_propagation',6*H*Hdot-Q*g*Qdot)
    zero('contracting_clock_response',(-Qdot/Q**2)/Hdot+6*H/(Q**3*g))
    # Sourced stationary hierarchy, evaluated exactly at the declared scale.
    ck=s.Integer(10)**30;de=s.Rational(1,10**140);LH=s.Rational(1,10**6)
    ell=1/ck;eta=1/(2*ck);tau=de**2/ck**2;acc=s.Rational(1,100);gam=s.Rational(11,10)
    mu=1-1/((1+acc)*(1+acc**2/100))
    r2=(ell+8*eta)/LH**2;r4=tau/LH**4
    relative=gam*((2*ell+8*eta)/LH**2+((ell+8*eta)**2+3*tau)/LH**4)/mu
    checks['stationary_scale_margin']=bool(relative<s.Rational(1,10**14))
    assert checks['stationary_scale_margin']
    # Ideal-fluid cubic density in independent current variables, with N_i=partial_i psi.
    e,aa,zz,nn,ww,asq,j2,jb,b2=s.symbols('e alpha zeta d w a2 j2 jb b2')
    vel=j2+2*s.exp(-2*e*zz)*(1+e*nn)*jb+s.exp(-4*e*zz)*(1+e*nn)**2*b2
    core=(1+e*nn)**2-e**2*asq*s.exp(2*e*zz)/(1+e*aa)**2*vel
    full=(1+e*aa)*s.exp(-3*ww*e*zz)*core**((1+ww)/2)
    exact3=s.diff(full,e,3).subs(e,0)/6
    p3=(-s.Rational(9,2)*ww**3*zz**3+s.Rational(9,2)*ww**2*aa*zz**2
        +(1+ww)*nn*(s.Rational(9,2)*ww**2*zz**2-3*ww*aa*zz)
        +ww*(1+ww)*nn**2*(aa-3*ww*zz)/2+ww*(1+ww)*(ww-1)*nn**3/6)
    proposed=p3-(1+ww)*asq/2*(((ww-1)*nn+(2-3*ww)*zz-aa)*(j2+2*jb+b2)+2*(nn-2*zz)*(jb+b2))
    zero('complete_ideal_fluid_cubic_density',exact3-proposed)
    report={'all_passed':all(checks.values()),'checks':checks,
      'endpoint':{'K':str(Kd),'V_numerator':str(s.expand(nV)),'V_denominator':str(denV),
                  'dV_numerator':str(s.expand(nv)),'dV_denominator':str(dv),
                  'kinetic_dissipation_numerator':str(s.expand(nk)),'kinetic_dissipation_denominator':str(dk)},
      'stationary':{'Hd_L':str(LH),'s':str(acc),'mu':str(mu),'second_derivative_ratio':str(r2),'sixth_derivative_ratio':str(r4),'relative_correction_bound':str(relative),
        'relative_correction_certified_upper':'1/10^14'},
      'limitations':['Does not prove nonlinear PDE well posedness','Does not compute a complete quantum cutoff','Does not construct a galaxy matched to the cosmological clock']}
    (ROOT/'data/exact_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'{len(checks)} symbolic and exact scale checks passed')

if __name__=='__main__':main()
