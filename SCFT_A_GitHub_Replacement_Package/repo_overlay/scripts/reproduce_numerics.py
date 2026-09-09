#!/usr/bin/env python3
"""Recompute the revised coefficient point; never reuse old response arrays.

CSV decimals are numerical approximations. Exact parameters and theorem
checks live in coefficient_certificate.json and exact_checks.json.
The complement of the flat switch is evaluated directly: subtracting its
value from one would lose the terms multiplied by kappa_m/Hd^2 = 10^30.
"""
import argparse, csv, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import expit
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'
KAP=1e30; DELTA=1e-140; QD=1.03; BD=80/103; KD=-BD*BD/10
QR=1.93; QE=1.99; QP=1.0375; WIDTH=QR-QP
ELL=1/(2*KAP); ETA=ELL; NU=DELTA/(2*KAP)
# The exact value is 10^-340. Binary64 would underflow it to zero;
# retain its exponent in long double before assembling the ODE updates.
TAU=np.longdouble('1e-340')
assert TAU>0, 'An extended exponent range is required for tau_6'
RHO0=np.array([.15,.000387,.000258]); WFL=np.array([0.,1/3,1/3])

def step(z):
    if z<=0:return 0.
    if z>=1:return 1.
    return float(expit(-1/z+1/(1-z)))

def rho_tail(q):
    v=step((q-QR)/(QE-QR))
    return (1-v)/((q-1)*(1-(q-1)**2))+8*v

def tail_solution():
    r=QR-1
    ini=[r/math.sqrt(1-r*r),r*r/(1+math.sqrt(1-r*r))]
    return solve_ivp(lambda q,y:[rho_tail(q)*y[0],y[0]],(QR,QE),ini,
                     method='DOP853',rtol=2e-13,atol=2e-15,dense_output=True,max_step=.001)

TAIL=tail_solution()
FE,KE=TAIL.y[:,-1]

def functions(q):
    """Return K,Kq,Kqq,Bq,Bqq,Sigma, in Hd units; q >= Qd."""
    if q>=QR:
        if q<QE:
            f,k=TAIL.sol(q);g=rho_tail(q)*f
        else:
            ex=math.exp(8*(q-QE));f=FE*ex;k=KE+FE*math.expm1(8*(q-QE))/8;g=8*f
        return KAP*k,KAP*f,KAP*g,0.,0.,1/16
    r=q-1;root=math.sqrt(1-r*r)
    db=KAP*r*r/(1+root);db1=KAP*r/root;db2=KAP/root**3
    u=q-QD
    kd=-3-3*BD*u+KD*u*u/2;kd1=-3*BD+KD*u
    if q<=QP:
        comp=cp=cpp=0.
    else:
        z=(QR-q)/WIDTH
        comp=float(expit(1/z-1/(1-z)))
        gp=(1/z**2+1/(1-z)**2)/WIDTH
        gpp=(2/z**3-2/(1-z)**3)/WIDTH**2
        cp=comp*(1-comp)*gp
        cpp=comp*(1-comp)*((1-2*comp)*gp*gp+gpp)
    k=kd+(db-kd)*comp
    k1=kd1+(db1-kd1)*comp+(db-kd)*cp
    k2=KD+(db2-KD)*comp+2*(db1-kd1)*cp+(db-kd)*cpp
    b1=BD*((1-comp)-u*cp);b2=-BD*(2*cp+u*cpp)
    sigma=1-(15/16)*step((q-1.75)/.15)
    return k,k1,k2,b1,b2,sigma

def background_at(N,loaded=True):
    F=math.exp(-3*N);rho=RHO0*np.exp(-3*(1+WFL)*N) if loaded else np.zeros(3)
    rv=float(sum(rho));hv=float(sum((1+WFL)*rho))
    def get(q):
        k,k1,k2,b1,b2,sigma=functions(q)
        rad=(q*b1)**2+4*(rv+q*k1-k)/3
        if rad<=0:raise ValueError(('negative Friedmann discriminant',N,q,rad))
        H=(q*b1+math.sqrt(rad))/2
        return k,k1,k2,b1,b2,sigma,H
    hi=1.2
    while get(hi)[1]+3*get(hi)[6]*get(hi)[3] < F:hi+=.3
    def current(q):
        k,k1,k2,b1,b2,sigma,H=get(q)
        return (k1+3*H*b1-F)/max(1.,F)
    q=brentq(current,QD,hi,xtol=5e-16,rtol=9e-16)
    k,k1,k2,b1,b2,sigma,H=get(q)
    ds=k2+3*H*b2;xi=ds+1.5*b1*b1
    qp=(-6*H*F+3*b1*(hv+q*F))/(2*H*xi)
    hp=-(ds*(hv+q*F)+3*H*b1*F)/(2*H*xi)
    fres=(3*H*H-rv-q*F+k)/max(1.,3*H*H,rv,abs(q*F),abs(k))
    cres=(k1+3*H*b1-F)/max(1.,abs(k1),abs(3*H*b1),F)
    return np.array([N,math.exp(N),q,H,qp,hp,ds,xi,b1,sigma,fres,cres])

COLS=['ln_a','a','Q','H_over_Hd','Q_prime','H_prime_over_Hd','D_over_Hd2','Xi_over_Hd2','Bq_over_Hd','Sigma_A','friedmann_scaled_residual','current_scaled_residual']
def csvwrite(path,header,rows):
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(header);w.writerows(rows)

def make_background(nodes):
    ns=np.linspace(0,math.log(10),nodes)
    bg=np.array([background_at(n) for n in ns])
    return bg,CubicSpline(ns,bg[:,2:10],axis=0)

def make_rhs(sp,k):
    def coefficients(N):
        Q,H,Qp,Hp,Ds,Xi,Bq,sigma=sp(N)
        x=k*k*math.exp(-2*N);r2=x/(H*H);bh=NU*x;A=2+3*bh
        j=2-Q*Bq/H+3*bh
        den=2*Q*Q*Xi/(H*H)+4*sigma*r2+3*bh*(Q*Q*Ds/(H*H)+2*sigma*r2)
        rhob=RHO0*np.exp(-3*(1+WFL)*N)/(H*H);hb=(1+WFL)*rhob
        return x,r2,bh,A,j,den,rhob,hb,Hp/H
    def rhs(N,y):
        x,r2,b,A,j,den,rho,hb,eh=coefficients(N)
        z,u=y[:2];dens=y[2:5];vel=y[5:8]
        mv=hb@vel;md=rho@dens
        alpha=(2*j*r2*u+3*j*mv+A*md-2*A*r2*(1+ELL*x)*z)/den
        zp=(b*r2*u+j*alpha-mv)/A
        up=(eh-1)*u+(1+ELL*x)*alpha+(1-8*ETA*x-3*TAU*x*x)*z
        dp=-(1+WFL)*(3*zp+r2*(vel-u))
        vp=(eh+3*WFL)*vel+alpha+WFL/(1+WFL)*dens
        return np.asarray(np.r_[zp,up,dp,vp],dtype=float)
    return rhs,coefficients

def channels(N,y,rhs,coeff):
    arr=[]
    for n,v in zip(N,y.T):
        x,r2,b,A,j,den,rho,hb,eh=coeff(n)
        z,u=v[:2];md=rho@v[2:5];mv=hb@v[5:8]
        al=(2*j*r2*u+3*j*mv+A*md-2*A*r2*(1+ELL*x)*z)/den
        yp=rhs(n,v)
        psi=u-z;phi=al-yp[1]+eh*u;dn=v[2:5]+3*(1+WFL)*u
        slip=phi-psi+ELL*x*al-(8*ETA*x+3*TAU*x*x)*z
        arr.append(np.r_[v,al,phi,psi,dn,slip,den])
    return np.array(arr,dtype=float)

def response(sp,k,fine=False,method='DOP853'):
    rhs,coeff=make_rhs(sp,k)
    ns=np.linspace(0,math.log(10),1001);ic=np.array([0,0,1,0,0,0,0,0.])
    sol=solve_ivp(rhs,(ns[0],ns[-1]),ic,method=method,t_eval=ns,
                  rtol=2e-11 if fine else 2e-10,atol=2e-16 if fine else 2e-15,
                  max_step=.005 if fine else .01)
    if not sol.success:raise RuntimeError(sol.message)
    return ns,channels(ns,sol.y,rhs,coeff)

def endpoint_and_contraction():
    mp.mp.dps=360 # resolves delta corrections even when the leading term is finite.
    d=mp.mpf(1)/10**140;c=mp.mpf(10)**30
    def endpoint(x):
        zn=180*c+125*d*x*x+371*d*x
        K=(500*c*x+224*c+375*d*x*x-12*d*x)/zn
        B=lambda y:2*(mp.mpf(6)/5+3*d*y/(2*c))*y*(1+y/(2*c))/( (mp.mpf(6)/5)**2+d*y/(2*c)*(2*y+mp.mpf(742)/125))
        D=(mp.mpf(6)/5)**2+d*x/(2*c)*(2*x+mp.mpf(742)/125)
        V=d*x*x*x*(1+x/(2*c))**2/(c*D)-x+4*x*x/c+3*d*d*x**3/c**2+3*B(x)-2*x*mp.diff(B,x)
        return K,V
    rows=[]
    for i in range(301):
        x=mp.mpf(1)/1000+(mp.mpf(250)-mp.mpf(1)/1000)*i/300
        K,V=endpoint(x)
        rows.append([i,mp.nstr(x,85),mp.nstr(K,85),mp.nstr(V,85),mp.nstr(V-(mp.mpf(2)/3*x+mp.mpf(19)/6*x*x/c),85)])
    csvwrite(DATA/'endpoint_spectrum.csv',['index','X_approx','K_approx','V_over_Hd2_approx','V_minus_delta_zero_limit_approx'],rows)
    # Charge is initial data: choose C_T=kappa_m for a transparent finite-scale-factor contracting example.
    rows=[]
    for i in range(101):
        q=mp.mpf(101)/100+mp.mpf(1)/200*i/100;r=q-1
        F=r/mp.sqrt(1-r*r);g=(1-r*r)**(-mp.mpf(3)/2);u=1-mp.sqrt(1-r*r)
        aa=F**(-mp.mpf(1)/3);hh=mp.sqrt((q*F-u)/3)
        qdot=3*hh*F/g
        rows.append([i,f'101/100 + {i}/20000',mp.nstr(q,85),mp.nstr(aa,85),mp.nstr(-hh,85),mp.nstr(qdot,85),mp.nstr(-6*(-hh)/(q**3*g),85)])
    csvwrite(DATA/'contracting_branch.csv',['index','Q_exact','Q_approx','a_approx','H_over_sqrt_kappa_approx','Qdot_over_sqrt_kappa_approx','sqrt_kappa_dNT_dH_approx'],rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--background-only',action='store_true');args=ap.parse_args()
    DATA.mkdir(exist_ok=True)
    broad=np.linspace(-15*math.log(10),math.log(900),601)
    for load,name in [(False,'source_free_background.csv'),(True,'visible_loaded_background.csv')]:
        bg=np.array([background_at(n,load) for n in broad]);csvwrite(DATA/name,COLS,bg)
        print(name,'Q range',bg[:,2].min(),bg[:,2].max(),'Xi min',bg[:,7].min(),'Qprime max',bg[:,4].max(),flush=True)
    bgc,spc=make_background(4001);bgf,spf=make_background(8001)
    csvwrite(DATA/'loaded_background.csv',COLS,bgf)
    if args.background_only:return
    endpoint_and_contraction()
    ks=.03*(100000/3)**(np.arange(42)/41)
    series=[];summary=[]
    labels=['zeta','u','delta_b','delta_gamma','delta_nu','theta_b','theta_gamma','theta_nu','alpha','Phi_N','Psi_N','delta_b_N','delta_gamma_N','delta_nu_N','slip_residual','canonical_denominator']
    for j,k in enumerate(ks):
        ns,co=response(spc,k);_,fi=response(spf,k,True)
        err=float(np.max(np.max(abs(fi[:,:14]-co[:,:14]),axis=0)/np.maximum(1,np.max(abs(fi[:,:14]),axis=0))))
        mx=float(abs(fi[:,[0,9,10,11,12,13]]).max())
        summary.append([j,float(k),mx,err,float(abs(fi[:,14]).max()),float(fi[:,15].min())]);series.append(fi)
        if j in [0,24,30,35,38,41]:csvwrite(DATA/f'response_{j:02d}.csv',['ln_a',*labels],np.column_stack([ns,fi]))
        print('mode',j,'maximum',repr(mx),'refinement',repr(err),flush=True)
    independent=[]
    for k in [.03,10.,1000.]:
        _,one=response(spf,k,True);_,two=response(spf,k,True,'RK45')
        err=float(np.max(np.max(abs(one[:,:14]-two[:,:14]),axis=0)/np.maximum(1,np.max(abs(one[:,:14]),axis=0))))
        independent.append({'k_over_Hd':k,'channel_scaled_difference':err})
    np.savez_compressed(DATA/'all_responses.npz',ln_a=ns,k_over_Hd=ks,channels=np.array(labels),responses=np.array(series))
    csvwrite(DATA/'response_summary.csv',['mode','k_over_Hd_approx','maximum_sampled_unit_response','refinement_discrepancy','slip_residual','min_canonical_denominator'],summary)
    report={'coefficient_point':'SCFT-A scale-separated revision, 2026-09-08',
      'exact_parameters':{'kappa_m_over_Hd2':'10^30','delta':'1/10^140','sigma':'1/10^6','s_N':'10','a0_over_Hd':'1'},
      'max_sampled_response':max(r[2] for r in summary),'max_refinement_difference':max(r[3] for r in summary),
      'max_slip_residual':max(r[4] for r in summary),'min_canonical_denominator':min(r[5] for r in summary),
      'max_background_scaled_residual':float(abs(bgf[:,10:12]).max()),
      'independent_integrators':independent,'tail_FE_over_kappa':float(FE),'tail_KE_over_kappa':float(KE),
      'tau6_Hd4_exact':'1/10^340',
      'tau6_storage':'numpy.longdouble; positive, with extended exponent range',
      'precision_note':'ODE updates are binary64. For 1<=a<=10 and k/Hd<=1000, 3*tau6*p^4 <= 3e-328; its effects are below binary64 update resolution. Exact high-momentum statements are symbolic, not inferred from this scan.',
      'scope':'Finite-grid numerical approximations; exact proofs and interval bounds are in the other reports.'}
    (DATA/'numerical_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
