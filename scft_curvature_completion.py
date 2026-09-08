"""Exact endpoint certificate and reproducible corrected SCFT-A response.

Units: H_d=M_Pl=1 in the numerical equations only.  The action coefficient
eta=1/(2*kappa_m) is exact.  delta_num is the uncertified numerical
evaluation printed in the supplied source; endpoint positivity is proved
for every 0<delta<10^-6 and does not depend on this decimal evaluation.
"""
from pathlib import Path
import json, argparse
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit
from scipy.interpolate import CubicSpline
import sympy as s

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'results')
OUT=parser.parse_args().output.resolve()
OUT.mkdir(parents=True,exist_ok=True)
checks={}
def check(name, value):
    checks[name]=bool(value)
    if not value:
        raise AssertionError(name)

X,d=s.symbols('X delta', positive=True)
Z=125*d*X**2+371*d*X+21600
D=Z/15000
B=2*(s.Rational(6,5)+d*X/80)*X*(1+X/240)/D
K=3*(125*d*X**2+(20000-4*d)*X+8960)/Z
V=(d*X/120)*X**2*(1+X/240)**2/D-X+d**2*X**3/4800+3*B-2*X*s.diff(B,X)
Vp=s.factor(V+X**2/30)
poly=s.Poly(s.cancel(57600*Z**2*Vp/X),X)
check('endpoint_positive_polynomial_identity',s.cancel(Vp-X*poly.as_expr()/(57600*Z**2))==0)
check('all_polynomial_coefficients_positive',all(all(a>0 for a in s.Poly(c,d).all_coeffs() if a!=0) for c in poly.all_coeffs()))
check('endpoint_auxiliary_denominator_positive',all(a>0 for a in s.Poly(Z,X).all_coeffs()))
check('endpoint_kinetic_positive_on_delta_interval',(20000-4*s.Rational(1,10**6))>0)
check('fixed_momentum_small_delta_limit',s.simplify(s.limit(Vp,d,0)-X*(19*X+480)/720)==0)
check('old_witness_limit',s.limit(V.subs(X,120),d,0)==-20)
check('corrected_witness_limit',s.limit(Vp.subs(X,120),d,0)==460)
witness=s.factor(Vp.subs(X,120))
wn,wd=s.fraction(s.factor(witness-460))
check('corrected_witness_strictly_above_460',all(c>0 for c in s.Poly(wn,d).all_coeffs() if c!=0) and all(c>0 for c in s.Poly(wd,d).all_coeffs() if c!=0))

# Independent auxiliary elimination of the complete scalar quadratic form.
ak,ws,h,bb,cc,zz,al,chi,zeta,eta,x,tau=s.symbols('ak ws h b C z alpha chi zeta eta x tau')
lag=-3*zz**2+ws*al*(3*zz-chi)+2*zz*chi+ak*al**2/2+2*cc*al*zeta+(x-8*eta*x**2-3*tau*x**3)*zeta**2-bb*(3*zz-chi-3*h*al)**2/2
DD=ws**2+bb*(ak+6*h*ws)
JJ=ws+3*bb*h
KK=((2+3*bb)*ak+3*ws**2+18*bb*h*(ws-h))/DD
solal=(2*JJ*zz-2*bb*cc*zeta)/DD
solchi=KK*zz+2*JJ*cc*zeta/DD
check('full_lapse_constraint',s.factor(s.diff(lag,al).subs({al:solal,chi:solchi}))==0)
check('full_shift_constraint',s.factor(s.diff(lag,chi).subs({al:solal,chi:solchi}))==0)
lr=s.factor(lag.subs({al:solal,chi:solchi}))
check('curvature_term_leaves_kinetic_coefficient',s.factor(s.diff(lr,zz,2)/2-KK)==0)
check('curvature_term_restoring_shift',s.factor(s.diff(lr,zeta,2)/2-(x-8*eta*x**2-3*tau*x**3-2*bb*cc**2/DD))==0)
check('auxiliary_hessian_determinant',s.factor(s.det(s.hessian(lag,(al,chi)))+DD)==0)
n,R,R1,R2,aa,ze,eps=s.symbols('N R R1 R2 alpha zeta eps')
check('new_hamiltonian_lapse_hessian_zero',s.diff(n*eta*R**2/2,n,2)==0)
cubic=s.expand(-eta*s.exp(3*eps*ze).series(eps,0,2).removeO()*(1+eps*aa)*(eps*R1+eps**2*R2)**2/2).coeff(eps,3)
check('intrinsic_curvature_cubic_vertex',s.expand(cubic+eta*(2*R1*R2+(aa+3*ze)*R1**2)/2)==0)
check('new_operator_no_velocity',s.diff(-n*eta*R**2/2,zz)==0)
report={'checks':checks,'number_passed':sum(checks.values()),'eta_Hd_squared':'1/240','endpoint_polynomial':str(poly.as_expr()),'endpoint_coefficients_ascending':[str(poly.nth(j)) for j in range(7)],'witness_exact':str(witness),'witness_margin_above_460':str(s.factor(witness-460))}
(OUT/'exact_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print('Exact checks passed:',len(checks),flush=True)

delta=float('8.0835430681749942437868886561856776507281743783432e-60')
km=120.; qd=1.03; bd=.8/qd; kd=-bd**2/10
nu=delta/240.; ell=1/240.; tau6=delta**2/14400.; eta0=1/240.
rh0=np.array([.15,.000387,.000258]); w=np.array([0.,1/3,1/3])

def switch(t):
    if t<=0:return 0.,0.,0.
    if t>=1:return 1.,0.,0.
    ss=expit(-1/t+1/(1-t))
    l=1/t**2+1/(1-t)**2
    return ss,ss*(1-ss)*l,ss*(1-ss)*((1-2*ss)*l*l-2/t**3+2/(1-t)**3)

def functions(q):
    # All numerical response backgrounds remain above qd: lower window=1.
    assert q>=qd-1e-12 and q<1.93
    v,v1,v2=switch((1.93-q)/.8925)
    v1/=-.8925;v2/=.8925**2
    r=q-1.; rt=np.sqrt(1-r*r)
    db=km*r*r/(1+rt); db1=km*r/rt; db2=km/rt**3
    dq=q-qd
    kk=-3-3*bd*dq+kd*dq*dq/2;kk1=-3*bd+kd*dq
    kval=(1-v)*db+v*kk
    k1=(1-v)*db1+v*kk1+v1*(kk-db)
    k2=(1-v)*db2+v*kd+2*v1*(kk1-db1)+v2*(kk-db)
    return kval,k1,k2,bd*(v+dq*v1),bd*(2*v1+dq*v2)

def background(N):
    f=np.exp(-3*N);rho=rh0*np.exp(-3*(1+w)*N)
    def equation(q):
        kv,k1,k2,b1,b2=functions(q)
        HH=np.sqrt((rho.sum()+q*f-kv)/3)
        return k1+3*HH*b1-f
    q=brentq(equation,qd,1.25,xtol=5e-15,rtol=1e-14)
    kv,k1,k2,b1,b2=functions(q)
    HH=np.sqrt((rho.sum()+q*f-kv)/3)
    ds=k2+3*HH*b2;Xi=ds+1.5*b1*b1
    hv=np.dot(1+w,rho)
    qp=(-6*HH*f+3*b1*(hv+q*f))/(2*HH*Xi)
    hp=-(ds*(hv+q*f)+3*HH*b1*f)/(2*HH*Xi)
    residual=max(abs(equation(q))/max(1,abs(k1),abs(3*HH*b1),f),abs(3*HH*HH-rho.sum()-q*f+kv)/max(1,3*HH*HH))
    return np.array([q,HH,qp,hp,ds,Xi,2*HH-q*b1,residual])

def build_background(nodes):
    grid=np.linspace(0,np.log(10),nodes)
    vals=np.array([background(N) for N in grid])
    return grid,vals,CubicSpline(grid,vals[:,:7],axis=0)

coarsegrid,coarsevals,bc=build_background(4001)
finegrid,finevals,bf=build_background(8001)
np.savetxt(OUT/'loaded_background.csv',np.column_stack((finegrid,finevals)),delimiter=',',header='ln_a,Q,H,Q_prime,H_prime,D_s,Xi_s,W_s,scaled_constraint_residual',comments='')
grid=np.linspace(0,np.log(10),1001)
initial=np.array([0.,0.,1.,0.,0.,0.,0.,0.])

def coefficients(N,k,bg):
    q,HH,qp,hp,ds,Xi,ws=bg(N)
    x=k*k*np.exp(-2*N);r2=x/(HH*HH);b=nu*x;A=2+3*b
    J=ws/HH+3*b
    T=2*q*q*Xi/(HH*HH)+4*r2+3*b*(q*q*ds/(HH*HH)+2*r2)
    rho=rh0*np.exp(-3*(1+w)*N)/(HH*HH)
    return HH,hp/HH,x,r2,b,A,J,T,rho

def rhs(N,y,k,bg):
    HH,eh,x,r2,b,A,J,T,rho=coefficients(N,k,bg)
    z,u=y[:2];dens=y[2:5];vel=y[5:8]
    mom=np.dot((1+w)*rho,vel);dd=np.dot(rho,dens)
    alpha=(2*J*r2*u+3*J*mom+A*dd-2*A*r2*(1+ell*x)*z)/T
    zp=(b*r2*u+J*alpha-mom)/A
    up=(eh-1)*u+(1+ell*x)*alpha+(1-8*eta0*x-3*tau6*x*x)*z
    denp=-(1+w)*(3*zp+r2*(vel-u))
    velp=(eh+3*w)*vel+alpha+w/(1+w)*dens
    return np.r_[zp,up,denp,velp]

def integrate(k,bg,method='DOP853',fine=False):
    sol=solve_ivp(lambda N,y:rhs(N,y,k,bg),(0,np.log(10)),initial,method=method,t_eval=grid,rtol=2e-11 if fine else 2e-10,atol=2e-16 if fine else 2e-15,max_step=.005 if fine else .01)
    if not sol.success:raise RuntimeError(sol.message)
    vals=[];res=[]
    for N,y in zip(grid,sol.y.T):
        HH,eh,x,r2,b,A,J,T,rho=coefficients(N,k,bg)
        z,u=y[:2];dens=y[2:5];vel=y[5:8]
        mom=np.dot((1+w)*rho,vel);dd=np.dot(rho,dens)
        alpha=(2*J*r2*u+3*J*mom+A*dd-2*A*r2*(1+ell*x)*z)/T
        dy=rhs(N,y,k,bg)
        psi=u-z;phi=alpha-dy[1]+eh*u
        dn=dens+3*(1+w)*u
        vals.append(np.r_[z,phi,psi,dn,y,alpha])
        q,HH,qp,hp,ds,Xi,ws=bg(N)
        ak=q*q*ds/(HH*HH)-6+6*(2-ws/HH)+2*r2
        E=ak-9*b
        terms=[E*alpha,-J*r2*u,3*J*dy[0],2*r2*(1+ell*x)*z,-dd]
        shift=[b*r2*u,-A*dy[0],J*alpha,-mom]
        slip=phi-psi+ell*x*alpha-(8*eta0*x+3*tau6*x*x)*z
        res.append([abs(sum(terms))/max(1,*map(abs,terms)),abs(sum(shift))/max(1,*map(abs,shift)),abs(slip)/max(1,abs(phi),abs(psi),abs(ell*x*alpha),abs((8*eta0*x+3*tau6*x*x)*z))])
    return np.array(vals),np.array(res),sol.nfev

ks=.03*(100000/3)**(np.arange(42)/41)
rows=[];data=[];norms=[]
for j,k in enumerate(ks):
    vc,rc,nc=integrate(k,bc)
    vf,rf,nf=integrate(k,bf,fine=True)
    ch=float(np.max(np.max(abs(vf-vc),axis=0)/np.maximum(1,np.max(abs(vf),axis=0))))
    pt=float(np.max(abs(vf-vc)/np.maximum(1,abs(vf))))
    obsmax=float(np.max(abs(vf[:,:6])))
    im=np.unravel_index(np.argmax(abs(vf[:,:6])),vf[:,:6].shape)
    rows.append([j,k,obsmax,np.exp(grid[im[0]]),im[1],ch,pt,*rf.max(axis=0),nf])
    data.append(vf)
    if j in (0,24,30,35,38,41):
        np.savetxt(OUT/f'response_{j:02}.csv',np.column_stack((grid,vf)),delimiter=',',header='ln_a,zeta_obs,Phi_N,Psi_N,delta_b_N,delta_gamma_N,delta_nu_N,zeta,u,delta_b,delta_gamma,delta_nu,vartheta_b,vartheta_gamma,vartheta_nu,alpha',comments='')
        print('mode',j,'maximum',repr(obsmax),'refinement',repr(ch),flush=True)
rows=np.array(rows)
np.savetxt(OUT/'response_summary.csv',rows,delimiter=',',header='j,k_Hd,maximum_observable_response,a_at_max,channel_at_max,channel_refinement,pointwise_refinement,lapse_residual,shift_residual,slip_residual,nfev',comments='')
np.savez_compressed(OUT/'all_responses.npz',N=grid,k=ks,values=np.array(data))
independent=[]
for k in (.03,10.,1000.):
    a,_,_=integrate(k,bf,fine=True)
    b,_,_=integrate(k,bf,method='RK45',fine=True)
    independent.append({'k_Hd':k,'channel_discrepancy':float(np.max(np.max(abs(a-b),axis=0)/np.maximum(1,np.max(abs(a),axis=0))))})
summary={'delta_num':str(delta),'delta_evaluation':'Inherited uncertified numerical evaluation from supplied TeX; no new extrema certificate. Exact endpoint proof covers the entire allowed interval.','eta_Hd_squared':'1/240','number_modes':42,'max_observable_response':float(rows[:,2].max()),'max_channel_refinement':float(rows[:,5].max()),'max_pointwise_refinement':float(rows[:,6].max()),'max_residuals':rows[:,7:10].max(axis=0).tolist(),'background_residual':float(finevals[:,-1].max()),'independent_integrators':independent}
(OUT/'numerical_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
# Plot values are evaluated from exact rational formulas converted to binary64.
ff=s.lambdify((X,d),Vp,'numpy');kk=s.lambdify((X,d),K,'numpy')
xx=np.unique(np.r_[np.linspace(.001,250,301),120.])
np.savetxt(OUT/'endpoint_spectrum.csv',np.column_stack((xx,kk(xx,delta),ff(xx,delta))),delimiter=',',header='p_squared_Hd_squared,K_d,V_d_Hd_squared',comments='')
print(json.dumps(summary,indent=2),flush=True)
