#!/usr/bin/python
import numpy as np
from scipy.optimize import *
import scipy as sp
import scipy.constants
import matplotlib.pyplot as plt
from matplotlib import rc
from sig_fig import *

class Parameter:
    def __init__(self, initialvalue, name):
        self.value = initialvalue
        self.name=name
    def set(self, value):
        self.value = value
    def __call__(self):
        return self.value

def fit(function, parameters, x, data, dx, dy):
    def fitfun(params):
        for i,p in enumerate(parameters):
            p.set(params[i])
        return (data - function(x))/np.sqrt(dy**2+(dx*params[0])**2) #DANGEROUS - only for use with affine fits
        
    if x is None: x = arange(data.shape[0])
    if dy is None: dy = ones(data.shape[0],"float")
    p = [param() for param in parameters]
    return leastsq(fitfun, p, full_output=1)

#Classical moment section
def clas_mom(plot=False, genTex = False):
    m_pos_tmp, curr, dcurr_tmp = np.loadtxt('classical_moment',unpack=True)
    dcurr = np.ones(curr.shape)*0.02 #correct the lazy entry of the data file
    m_pos = m_pos_tmp*0.1
    dm_pos = np.ones(m_pos.shape)*0.001
    
    b_field = 13.6*curr #in Gauss
    db_field = b_field*np.sqrt((0.3/13.6)**2+(dcurr/curr)**2)
    
    rmg = m_pos*1.39*980.438 #value of g from wolfram alpha
    drmg = rmg*np.sqrt((dm_pos/m_pos)**2+(0.01/1.39)**2)
    
    mu_inv = Parameter(1.,'mu_inv')
    yint = Parameter(1.,'yint')
    
    ps = [mu_inv,yint]
    p1,cov1,info1,msg1,success1 = fit((lambda x:x*mu_inv()+yint()),ps,rmg,b_field,drmg,db_field)
    chisq1 = sum(info1['fvec']*info1['fvec'])
    dof1 = len(b_field)-len(ps)
    
    mu = 1./mu_inv()
    dmu = mu*(np.sqrt(cov1[0,0])*np.sqrt(chisq1/dof1)/mu_inv())
    
    print 'slope=%s, yint=%s'%(sig_fig(mu_inv(),np.sqrt(cov1[0,0])*np.sqrt(chisq1/dof1)),sig_fig(yint(),np.sqrt(cov1[1,1])*np.sqrt(chisq1/dof1)))
    print 'mu = %s'%sig_fig(mu,dmu)
    print 'reduced chisq : %.3f'%(chisq1/dof1)

    if genTex:
        for i in range(len(curr)):
            print '$%.2f$ & $%s$\\\\'%(m_pos_tmp[i],sig_fig(curr[i],dcurr[i]))
    
    if plot:
        plt.figure()
        plt.errorbar(rmg,b_field,xerr=drmg,yerr=db_field,fmt='r.')
        indep1 = np.arange(rmg[0],rmg[len(rmg)-1:],300)
        plt.plot(indep1,indep1*mu_inv()+yint(),'b-')
        plt.gca().set_xlabel(r'rmg ($g\cdot\frac{\mathrm{cm}^2}{\mathrm{s}^2}$)')
        plt.gca().set_ylabel(r'Magnetic field (Gauss)')
        plt.savefig('grav-fit.png')
        plt.show()
    return [mu,dmu]

#Classical precession
def clas_pre(exp_mu,plot=False,genTex=False):
    curr, pd = np.loadtxt('classical_precession',unpack=True)
    dpd = 0.1*np.ones(pd.shape)
    dcurr = 0.02*np.ones(curr.shape)

    b_field = 13.6*curr #in Gauss
    db_field = b_field*np.sqrt((0.3/13.6)**2+(dcurr/curr)**2)


    omega = 2*sp.constants.pi/pd
    domega = omega*(dpd/pd)
    
    mul = Parameter(1., 'mul')
    yint = Parameter(1., 'yint')
    ps = [mul,yint]
    p1,cov,info,msg,success = fit((lambda x:x*mul()+yint()),ps,b_field,omega,db_field,domega)
    chisq = sum(info['fvec']*info['fvec'])
    dof = len(omega)-len(ps)
    dyint = np.sqrt(cov[1,1])*np.sqrt(chisq/dof)
    dmul = np.sqrt(cov[0,0])*np.sqrt(chisq/dof)

    print 'slope: %s\nyint: %s'%(sig_fig(mul(),dmul),sig_fig(yint(),dyint))
    
    mom_in = 2./5*137.5*(5.377/2)**2
    dmom_in = mom_in*np.sqrt((.01/137.5)**2+(.04/5.377)**2)
    print 'I=%s'%sig_fig(mom_in,dmom_in)

    ang_mom = mom_in*2*sp.constants.pi*5.1
    dang_mom = ang_mom*np.sqrt((dmom_in/mom_in)**2+(0.1/5.1)**2)
    
    print 'L=%s'%sig_fig(ang_mom,dang_mom)
    exp_coef = exp_mu[0]/ang_mom
    dexp_coef = exp_coef*np.sqrt((exp_mu[1]/exp_mu[0])**2+(dang_mom/ang_mom)**2)

    print 'Expected mu: %s'%sig_fig(exp_mu[0],exp_mu[1])
    print 'Got : %s with reduced chisq %.3f'%(sig_fig(mul()*ang_mom,mul()*ang_mom*np.sqrt((dmul/mul())**2+(dang_mom/ang_mom)**2)),chisq/dof)
    print 'y-int: %s'%sig_fig(yint(),dyint)
    print 'Residuals:'
    for i in range(len(pd)):
        print '%.3f: residual of %.3f'%(curr[i],omega[i]-b_field[i]*mul()+yint())
    print
    if genTex:
        for i in range(len(curr)):
            print '$%s$ & $%s$\\\\'%(curr[i],pd[i])
    if plot:
        plt.figure()
        plt.errorbar(b_field,omega,xerr=db_field, yerr=domega,fmt='g.')
        indep = np.arange(b_field[0],b_field[len(b_field)-1]+9,1)
        plt.plot(indep,indep*mul()+yint(),'k-')
        plt.gca().set_xlabel(r'B-field (Gauss)')
        plt.gca().set_ylabel(r'Precession frequency(rad/s)')
        plt.savefig('prec-fit.png')
        plt.show()

def electron_spin(plot=False, genTex=False):
    a = np.loadtxt('electron_flip')
    freqs = a[:,0]
    dfreqs = 0.01*np.ones(freqs.shape)

    currs = np.mean(a[:,1:],axis=1)
    sdcurrs = np.std(a[:,1:],axis=1)
    dcurrs = 0.04*np.ones(currs.shape)

    print sdcurrs
    bfield = currs*4.8 #gauss
    dbfield = 4.8*dcurrs

    ang_freqs = freqs*1000000*2*sp.constants.pi
    dang_freqs = dfreqs*1000000*2*sp.constants.pi

    gmr_inv = Parameter(1.,'gmr')
    yint = Parameter(1.,'yint')

    #ps = [gmr_inv,yint]
    ps = [gmr_inv]
    p1,cov,info,msg,success = fit((lambda x:x*gmr_inv()),ps,ang_freqs,bfield,dang_freqs,dbfield)
    chisq = sum(info['fvec']*info['fvec'])
    dof = len(bfield)-len(ps)
    #dyint = np.sqrt(cov[1,1])*np.sqrt(chisq/dof)
    gmr = 1/gmr_inv()
    dgmr = np.sqrt(cov[0,0])*np.sqrt(chisq/dof)/gmr_inv()**2

    print 'slope: %s'%sig_fig(gmr_inv(),np.sqrt(cov[0,0])*np.sqrt(chisq/dof))

    print 'GMR: %s'%sig_fig(gmr,dgmr)
    print gmr
    #print 'y int: %s'%sig_fig(yint(),dyint)
    print 'reduced chisq: %.3f'%(chisq/dof)

    if genTex:
        for i in range(len(freqs)):
            print '$%.2f$ & '%freqs[i],
            for volt in a[i,1:]:
                print '$%.2f$, '%volt,
            print '& $%.2f$\\\\'%currs[i]

    if plot:
        plt.figure()
        plt.errorbar(ang_freqs,bfield,xerr=dang_freqs,yerr=dbfield,fmt='g.')
        indep = np.arange(ang_freqs[0],ang_freqs[len(ang_freqs)-1:],1)
        plt.plot(indep,indep*gmr_inv(),'k-')
        plt.gca().set_xlabel('Resonance frequency(rad/s)')
        plt.gca().set_ylabel('Magnetic field(Gauss)')
        plt.savefig('elec-fit.png')
        plt.show()



#exp_mu = clas_mom(True)
#clas_pre(exp_mu,True,True)
electron_spin(True)
