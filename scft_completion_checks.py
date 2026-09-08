"""Exact general braided and matter auxiliary checks for the integrated SCFT action.

No nonlinear existence or no-loss energy theorem is inferred from these identities.
Run: python3 scft_completion_checks.py
"""
import json
from pathlib import Path
import sympy as s
b,H,W,Ak,C,x,tau,z,Z,al,chi,Vm,Dr=s.symbols('b H W Ak C x tau z Z alpha chi V_m Delta_rho', real=True)
A=2+3*b
J=W+3*b*H
D=W**2+b*(Ak+6*H*W)
K=(A*Ak+3*W**2+18*b*H*(W-H))/D
B=2*J*C/D
V0=2*b*C**2/D-x+3*tau*x**3
L=-3*z**2+W*al*(3*z-chi)+2*z*chi+Ak*al**2/2+2*C*al*Z+(x-3*tau*x**3)*Z**2-b*(3*z-chi-3*H*al)**2/2
sol={al:(2*J*z-2*b*C*Z)/D,chi:K*z+2*J*C*Z/D}
checks={}
def zero(name,e):
    r=s.factor(s.cancel(e))
    checks[name]={'passed':r==0,'residual':str(r)}
zero('general_lapse_auxiliary',s.diff(L,al).subs(sol))
zero('general_shift_auxiliary',s.diff(L,chi).subs(sol))
zero('general_auxiliary_determinant',s.det(s.hessian(L,(al,chi)))+D)
zero('general_reduced_lagrangian',L.subs(sol)-(K*z**2+2*B*z*Z-V0*Z**2))
Lm=L+Vm*chi-Dr*al-3*Vm*z
solm={al:(2*J*z-2*b*C*Z+J*Vm+b*Dr)/D,
      chi:K*z+2*J*C*Z/D+((Ak-9*b*H**2)*Vm-J*Dr)/D}
zero('matter_lapse_auxiliary',s.diff(Lm,al).subs(solm))
zero('matter_shift_auxiliary',s.diff(Lm,chi).subs(solm))
zero('base_lapse_limit',solm[al].subs({b:0,C:x})-(2*z+Vm)/W)
zero('base_shift_limit',solm[chi].subs({b:0,C:x})-((3+2*Ak/W**2)*z+2*x*Z/W+Ak*Vm/W**2-Dr/W))
M=s.symbols('M',real=True)
Du=4*H**2+b*(M+6*H**2+2*x)
zero('local_unbraided_D',D.subs({W:2*H,Ak:M-6*H**2+2*x})-Du)
zero('local_unbraided_K',K.subs({W:2*H,Ak:M-6*H**2+2*x})-A*(M+2*x)/Du)
nu,ell,Sigma=s.symbols('nu ell Sigma',positive=True)
zero('completed_general_UV_K',s.limit(K.subs({b:nu*x,Ak:M+2*Sigma*x}),x,s.oo)-3)
zero('completed_general_UV_V',s.limit(V0.subs({b:nu*x,Ak:M+2*Sigma*x,C:x*(1+ell*x)})/x**3,x,s.oo)-(ell**2/Sigma+3*tau))
# Exact stationary coupled determinant: positive scalar restoring, no physical pole.
stat=s.Matrix([[x,-x*(1+ell*x)],[-x*(1+ell*x),x-3*tau*x**3]])
zero('stationary_coupled_determinant',stat.det()+x*(2*ell*x**2+(ell**2+3*tau)*x**3))
# Constant-coefficient homogeneous Schur check including the braid.
N,P,B0,nu0,p2=s.symbols('N P B0 nu0 p2',positive=True)
Bfun=s.Function('B'); Ufun=s.Function('U'); q=1/N
# Homogeneous K inverse has A=1 for zero Fourier mode.
Kaux=s.symbols('Kaux')
ha=N*((s.Rational(2,3)*P-Bfun(q))*Kaux+Kaux**2/3-Ufun(q))
ksol=-P+s.Rational(3,2)*Bfun(q)
hred=s.simplify(ha.subs(Kaux,ksol))
neg=-N**2*s.diff(hred,N,2)
expected=q*(s.diff(Ufun(q),N,2)*0) # replaced explicitly below using dummy Q
Q=s.symbols('Q',positive=True)
expected=(Q*(s.diff(Ufun(Q),Q,2)+s.diff(Bfun(Q),Q,2)*(-P+s.Rational(3,2)*Bfun(Q))+s.Rational(3,2)*s.diff(Bfun(Q),Q)**2)).subs(Q,1/N)
zero('homogeneous_braid_Schur',neg-expected)
report={'scope':'Exact completed-action quadratic and auxiliary identities only; nonlinear analytic closure is conditional.',
        'floating_point_used':False,'checks':checks,'passed_count':sum(v['passed'] for v in checks.values()),'total_count':len(checks)}
report['all_passed']=report['passed_count']==report['total_count']
Path('scft_completion_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
assert report['all_passed']

