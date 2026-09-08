"""BASE MODEL S_0 background and finite-k audit; not completed-action mode predictions.

Dependencies: numpy, scipy. Numerical outputs are approximations, never exact
decimal constants. All dimensionless parameters below represent exact ratios.
Run in an empty directory; writes scft_background.csv and scft_modes.npz.
"""
import json, math
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.interpolate import CubicSpline
from scipy.special import expit

QD=103/100; QR=193/100; QE=199/100
BD=80/103; KD=-BD**2/10; KM=120.; LAM=8.
RB=3/20; RR=645/10**6

def switch(x):
    """Value and first four derivatives, with exact flat plateaux."""
    if x<=0: return np.array([0.,0.,0.,0.,0.])
    if x>=1: return np.array([1.,0.,0.,0.,0.])
    t=1/(1-x)-1/x
    if abs(t)>500: return np.array([float(t>0),0.,0.,0.,0.])
    s=expit(t); u=s*(1-s)
    t1=(1-x)**-2+x**-2
    t2=2*(1-x)**-3-2*x**-3
    t3=6*(1-x)**-4+6*x**-4
    t4=24*(1-x)**-5-24*x**-5
    return np.array([s,u*t1,u*((1-2*s)*t1*t1+t2),
      u*((1-6*s+6*s*s)*t1**3+3*(1-2*s)*t1*t2+t3),
      u*((1-14*s+36*s*s-24*s**3)*t1**4+
         6*(1-6*s+6*s*s)*t1*t1*t2+
         (1-2*s)*(3*t2*t2+4*t1*t3)+t4)])

def mul(a,b):
    return np.array([sum(math.comb(n,j)*a[j]*b[n-j]
                         for j in range(n+1)) for n in range(len(a))])

def core(q):
    x=q-1; d=1-x*x
    db=np.array([KM*x*x/(1+np.sqrt(d)), KM*x/np.sqrt(d),
        KM/d**1.5, 3*KM*x/d**2.5, 3*KM*(1+4*x*x)/d**3.5])
    dq=q-QD
    kd=np.array([-3-3*BD*dq+KD*dq*dq/2,-3*BD+KD*dq,KD,0.,0.])
    lo=switch((q-1-.03/2)/(.03/4))/np.array([(.03/4)**n for n in range(5)])
    scale=-(.03*(30-1/4))**-1
    up=switch((QR-q)/(.03*(30-1/4)))*np.array([scale**n for n in range(5)])
    win=mul(lo,up)
    kval=db+mul(win,kd-db)
    bval=mul(win,np.array([BD*dq,BD,0.,0.,0.]))
    return kval,bval

def rho(q):
    if q>=QE: return np.array([8.,0.,0.])
    x=q-1; d=x*(1-x*x); d1=1-3*x*x; d2=-6*x
    r=1/d; r1=-d1/d**2; r2=2*d1*d1/d**3-d2/d**2
    v=switch((q-QR)/(QE-QR))[:3]/np.array([1.,QE-QR,(QE-QR)**2])
    return np.array([r,r1,r2])+mul(v,np.array([8-r,-r1,-r2]))

KR, BR=core(QR); FR=KR[1]
def tail_rhs(q,z):
    return [rho(q)[0],np.exp(z[0])]
MATCH=solve_ivp(tail_rhs,(QR,QE),[np.log(FR),KR[0]],rtol=3e-13,atol=3e-13,
                dense_output=True,max_step=.0001)
FE=np.exp(MATCH.y[0,-1]); KE=MATCH.y[1,-1]

def functions(q):
    if q<QR: k,b=core(q)
    else:
        if q<QE:
            z=MATCH.sol(q); f=np.exp(z[0]); k0=z[1]
        else:
            f=FE*np.exp(8*(q-QE)); k0=KE+(f-FE)/8
        r,r1,r2=rho(q)
        k=np.array([k0,f,r*f,(r*r+r1)*f,(r**3+3*r*r1+r2)*f])
        b=np.zeros(5)
    s=switch((q-7/4)/(3/20))[:3]/np.array([1.,3/20,(3/20)**2])
    sig=np.array([1.,0.,0.])-15/16*s
    return k,b,sig

def background(n,loaded=True):
    a=np.exp(n); f=np.exp(-3*n)
    rb=RB*f if loaded else 0.; rr=RR*np.exp(-4*n) if loaded else 0.
    rv=rb+rr; hv=rb+4*rr/3
    if f>=FE: q=QE+np.log(f/FE)/8
    elif f>=FR:
        q=brentq(lambda q:MATCH.sol(q)[0]-np.log(f),QR,QE,xtol=5e-15,rtol=1e-15)
    else:
        def eq(q):
            k,b=core(q); h=np.sqrt(max((rv+q*f-k[0])/3,0))
            return k[1]+3*h*b[1]-f
        q=brentq(eq,QD,QR,xtol=5e-15,rtol=1e-15)
    k,b,sig=functions(q); h=np.sqrt((rv+q*f-k[0])/3)
    d=k[2]+3*h*b[2]; xi=d+1.5*b[1]**2; w=2*h-q*b[1]
    qp=(-6*h*f+3*b[1]*(hv+q*f))/(2*h*xi)
    hp=-(d*(hv+q*f)+3*h*b[1]*f)/(2*h*xi)
    cs=-2*h*hp+h*qp*(b[1]+q*b[2])+h*q*b[1]-.5*q*q*b[1]**2
    aa=q*q*d-6*h*h+6*h*q*b[1]
    ups=1-2*sig[0]*hv/(w*w)
    mix=np.sqrt(max(aa*hv/(w*w*ups),0))/h
    return np.array([q,h,qp,hp,xi,w,cs,sig[0],ups,mix,
        q*q*xi/h**2,w/h,aa/h**2,cs/h**2,rb/h**2,rr/h**2,
        k[0],k[1],k[2],k[3],k[4],sig[1],sig[2]])

def run(npoints=20001,ks=None,rtol=2e-9,atol=2e-13,max_step=.01):
    ns=np.linspace(np.log(1e-10),np.log(10),npoints)
    bg=np.array([background(n) for n in ns]); spl=CubicSpline(ns,bg,axis=0)
    np.savetxt('scft_background.csv',np.column_stack([ns,bg]),delimiter=',',
        header='N,Q,H,Qprime,Hprime,Xi,W,C,Sigma,Upsilon,kmix_over_aH,D0_over_H2,W_over_H,A_over_H2,C_over_H2,rhob_over_H2,rhor_over_H2,K,K1,K2,K3,K4,Sigma1,Sigma2')
    ks=np.array([.1,1.,10.,100.,1000.]) if ks is None else np.array(ks)
    modes={}; summary=[]; wa=np.array([0.,1/3])
    for kval in ks:
        def rhs(n,y):
            q,h,qp,hp,xi,w,c,sig,ups,mix,d0,wb,ab,cb,rb,rr,*_=spl(n)
            r2=(kval/(np.exp(n)*h))**2
            d=d0+2*sig*r2; ak=ab+2*sig*r2; qs=4*d/wb**2
            fs=4*cb/wb**2; hb=np.array([rb,4*rr/3]); rh=np.array([rb,rr])
            z,pi=y[:2]; de=y[2:4]; th=y[4:6]
            vv=hb@th; dd=rh@de
            zp=(pi-2*ak/wb**2*vv+2*dd/wb)/qs
            alpha=(2*zp+vv)/wb; chi=pi/2+2*r2*z/wb
            return np.r_[zp,-(3+hp/h)*pi-r2*(fs*z-2*vv/wb),
                  -(1+wa)*(3*zp+r2*th-chi),
                  (hp/h+3*wa)*th+alpha+wa/(1+wa)*de]
        r2=(kval/(np.exp(ns[0])*bg[0,1]))**2
        y0=np.array([1.,-4*r2/3,r2/6,2*r2/9,r2/36,r2/18])
        sol=solve_ivp(rhs,(ns[0],ns[-1]),y0,method='Radau',rtol=rtol,atol=atol,
                      max_step=max_step,dense_output=True)
        if not sol.success: raise RuntimeError(sol.message)
        no=np.linspace(ns[0],ns[-1],6001); yy=sol.sol(no); bb=spl(no)
        r2s=(kval/(np.exp(no)*bb[:,1]))**2
        shift=-yy[1]/(2*r2s)-2*yy[0]/bb[:,11]
        psi=-(yy[0]+shift); dbn=yy[2]-3*shift; drn=yy[3]-4*shift
        pk=r2s*bb[:,8]-bb[:,12]*(bb[:,14]+4*bb[:,15]/3)/bb[:,11]**2
        crossings=int(np.count_nonzero(pk[1:]*pk[:-1]<0))
        at1=sol.sol(0); b1=spl(0); r21=(kval/b1[1])**2
        sh1=-at1[1]/(2*r21)-2*at1[0]/b1[11]
        maxamp=float(np.max(np.abs(np.array([yy[0],psi,dbn,drn]))))
        summary.append({'k':float(kval),'max_r':float(np.sqrt(r2s.max())),
           'crossings':crossings,'db1':float(at1[2]-3*sh1),'db10':float(dbn[-1]),
           'psi1':float(-(at1[0]+sh1)),'maxamp':maxamp,'nfev':sol.nfev})
        modes[str(kval)]=np.column_stack([no,yy[0],psi,dbn,drn,pk,np.sqrt(r2s)])
        print(json.dumps(summary[-1]),flush=True)
    np.savez_compressed('scft_modes.npz',**modes)
    out={'FE':FE,'KE':KE,'background_samples':npoints,
         'minimum_Upsilon':float(bg[:,8].min()),'maximum_kmix_aH':float(bg[:,9].max()),
         'minimum_Xi':float(bg[:,4].min()),'minimum_W':float(bg[:,5].min()),
         'minimum_C':float(bg[:,6].min()),'maximum_Qprime':float(bg[:,2].max()),
         'maximum_Hprime':float(bg[:,3].max()),'modes':summary}
    Path('scft_audit.json').write_text(json.dumps(out,indent=2))
    return out

if __name__=='__main__':
    print(json.dumps(run(ks=np.unique(np.r_[np.geomspace(.03,1000,42),.1,1,10,100]))),flush=True)

