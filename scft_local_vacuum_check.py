"""HISTORICAL BASE MODEL checks; trace-repair substitutions are not current-action checks.

Run after compiling the manuscript, which writes scft_audit.py.
All numerical outputs are approximations; arithmetic checks use exact rationals.
This script does not claim a nonlinear, MOND-preserving vacuum repair.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss
from fractions import Fraction as F
from itertools import product


def symbolic_checks():
    # Exact rational substitutions check the separately derived closed forms.
    # These are arithmetic cross-checks, not a computer-algebra general proof.
    count=0
    for x,mu2,sig,lam,v,z in product(
        [F(1,9),F(1),F(16)], [F(1,4),F(1),F(60)],
        [F(1,16),F(1,2),F(1)], [F(-4),F(-1),F(-1,3),F(1,7),F(1),F(3)],
        [F(-2,3),F(3,5)], [F(-5,7),F(2,3)]):
        A=mu2+sig*x
        chi=(2+3*lam)*v/lam
        al=-x*z/A
        L=-3*v*v+2*v*chi+A*al*al+2*x*al*z+x*z*z-lam*(3*v-chi)**2/2
        kinetic=(2+3*lam)/lam
        potential=x*((1-sig)*x-mu2)/A
        assert L == kinetic*v*v-potential*z*z
        assert 2*v+lam*(3*v-chi) == 0
        assert A*al+x*z == 0
        count+=1
    for rootq,beta,rr,x in product([F(1,3),F(1),F(4,5)],
           [F(1,100),F(2,3)], [F(3,2),F(2),F(3)], [F(1,9),F(1),F(16)]):
        q=rootq**2; qdot=-rr*beta*q*rootq
        K=rr*beta*beta*q; bq=-beta/rootq
        assert q*rr*beta*beta-K == 0
        assert K-bq*qdot == 0
        W=beta*rootq; Wdot=-rr*beta*beta*q/2
        Xi=F(3,2)*beta*beta/q
        C=-Wdot-W*W/2
        Qs=4*(q*q*Xi+2*x)/W**2
        Fs=4*C/W**2
        assert Qs == 6+8*x/(beta*beta*q)
        assert Fs == 2*(rr-1)
        rho=F(3,7)
        zstat=rr*rho/(2*(rr-1)*x)
        assert Fs*x*zstat == -2*rho*Wdot/W**2
        psidot=(-2*zstat+rho/x)*(-Wdot/W**2)
        assert psidot == -zstat
        count+=1
    for h,w,wd,x in product([F(1),F(3,2)], [F(2,3),F(-3,2)],
                            [F(-3,5),F(4,3)], [F(1,9),F(1),F(16)]):
        pp,aa,zz=F(7,3),F(5,4),F(-1,3)
        cs=-wd+h*w-w*w/2
        phi=2*cs*zz/w**2+h*pp/(2*aa**3*x)+2*wd*zz/w**2
        psi=(2*h/w-1)*zz+h*pp/(2*aa**3*x)
        assert phi == psi
        count+=1
    return {'exact_rational_substitution_checks_passed': count,
            'general_proofs': 'analytic derivations in manuscript; no CAS proof claimed',
            'rolling_stationary_G_over_Gstar': 'r/(r-1)'}


def load_background(path):
    spec = importlib.util.spec_from_file_location('scft_base', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def packet_checks(base, nq=96, tol=1e-10):
    nodes, weights = leggauss(nq)
    kvals = 24 + 8*nodes
    u = (kvals-16)/16
    raw = 8*weights*kvals**2*np.exp(-1/(u*(1-u)))
    weights = raw/raw.sum()
    end = np.log(2.)
    times = np.linspace(0,end,301)

    def coeff(n):
        bg = base.background(n,loaded=False)
        q,h,qp,hp,xi,w,c,sig = bg[:8]
        k,b,ss = base.functions(q)
        xp = kvals**2*np.exp(-2*n)
        d = q*q*xi+2*sig*xp
        qs = 4*d/w**2
        fs = 4*c/w**2
        wp = 2*hp-qp*(b[1]+q*b[2])
        xip = qp*k[3]+3*hp*b[2]+3*h*b[3]*qp+3*b[1]*b[2]*qp
        dp = 2*q*qp*xi+q*q*xip+2*ss[1]*qp*xp-4*sig*xp
        logqsp = dp/d-2*wp/w
        return q,h,hp,xi,w,c,xp,d,qs,fs,logqsp,wp

    def canonical(n,y):
        q,h,hp,xi,w,c,xp,d,qs,fs,lq,wp=coeff(n)
        zz,pp = y[:nq],y[nq:]
        return np.r_[pp/(np.exp(3*n)*qs*h),-np.exp(3*n)*fs*xp*zz/h]

    def second_order(n,y):
        q,h,hp,xi,w,c,xp,d,qs,fs,lq,wp=coeff(n)
        zz,vv = y[:nq],y[nq:]
        return np.r_[vv,-(3+hp/h+lq)*vv-c*xp/(d*h*h)*zz]

    initial=np.r_[np.ones(nq),np.zeros(nq)]
    can=solve_ivp(canonical,(0,end),initial,method='DOP853',t_eval=times,
                  rtol=tol,atol=tol/100,max_step=.005)
    sec=solve_ivp(second_order,(0,end),initial,method='DOP853',t_eval=times,
                  rtol=tol,atol=tol/100,max_step=.005)
    assert can.success and sec.success
    err = float(np.max(np.abs(can.y[:nq]-sec.y[:nq])))
    radius=np.linspace(0,2,401)
    bessel=np.sinc(np.outer(kvals,radius)/np.pi)
    psis=[]; alphas=[]; zetas=[]; gs=[]; constraints=[]; marg=[]
    for i,n in enumerate(times):
        q,h,hp,xi,w,c,xp,d,qs,fs,lq,wp=coeff(n)
        zz,pp=can.y[:nq,i],can.y[nq:,i]
        zd=pp/(np.exp(3*n)*qs)
        alpha=2*zd/w
        shift=-pp/(2*np.exp(3*n)*xp)-2*zz/w
        newton=(2*h/w-1)*zz+h*pp/(2*np.exp(3*n)*xp)
        chi=-xp*shift
        ak=d-1.5*w*w
        lapse=ak*alpha+3*w*zd-w*chi+2*xp*zz
        scale=1+np.abs(ak*alpha)+np.abs(3*w*zd)+np.abs(w*chi)+np.abs(2*xp*zz)
        constraints.append(np.max(np.abs(lapse)/scale))
        psis.append((weights*newton)@bessel)
        alphas.append((weights*alpha)@bessel)
        zetas.append((weights*zz)@bessel)
        gs.append((weights*shift)@bessel)
        # A Fourier L1 bound controls the supremum at every spatial point.
        marg.append([q,h,xi,w,c,np.min(qs),fs,
                     np.sum(weights*np.abs(alpha)),np.sum(weights*np.abs(newton)),
                     np.min(.5*(pp*pp/(np.exp(3*n)*qs)+np.exp(3*n)*fs*xp*zz*zz))])
    marg=np.array(marg)
    assert np.all(marg[:,[0,1,2,3,4,5,6,9]]>0)
    out={'quadrature_nodes':nq,'k_over_Hd_interval':[16,32],
         'a_interval':[1,2], 'amplitude':1e-6,
         'mode_formulation_max_absolute_difference':err,
         'scaled_lapse_constraint_max_residual':float(max(constraints)),
         'minimum_Q':float(marg[:,0].min()),
         'minimum_Xi_over_Hd2':float(marg[:,2].min()),
         'minimum_W_over_Hd':float(marg[:,3].min()),
         'minimum_C_over_Hd2':float(marg[:,4].min()),
         'minimum_Qs_over_Mpl2':float(marg[:,5].min()),
         'minimum_Fs_over_Mpl2':float(marg[:,6].min()),
         'spatial_lapse_bound_max_over_sampled_times':float(marg[:,7].max()*1e-6),
         'spatial_Newtonian_bound_max_over_sampled_times':float(marg[:,8].max()*1e-6)}
    assert out['spatial_lapse_bound_max_over_sampled_times'] < min(.5,(out['minimum_Q']-1)/2)
    return out, {'times':times,'radius':radius,'Psi':np.array(psis),
                 'alpha':np.array(alphas),'zeta':np.array(zetas),'shift':np.array(gs)}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--background',type=Path,required=True)
    args=parser.parse_args()
    result=symbolic_checks()
    base=load_background(args.background)
    result['base_packet'], data=packet_checks(base,96,1e-10)
    result['refined_packet'], refined=packet_checks(base,192,1e-12)
    result['quadrature_refinement_max_absolute_packet_difference']=float(
        np.max(np.abs(data['Psi']-refined['Psi']))*1e-6)
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    main()
