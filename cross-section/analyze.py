#!/usr/bin/python
import numpy as np
from pylab import *
from scipy.optimize import leastsq

energy,thickness,count,time = loadtxt('data',unpack=True,usecols=[0,1,2,3])

toBeRemoved = []
currenergy = 0
bgrate = 0

for (i,t) in enumerate(thickness):
    if t == -1:
        currenergy = energy[i]
        bgrate = count[i]/time[i]
        toBeRemoved.append(i)
    else:
        if energy[i] == currenergy:
            count[i] = count[i]-bgrate*time[i]

energy=np.delete(energy,toBeRemoved)
thickness=np.delete(thickness,toBeRemoved)
count=np.delete(count,toBeRemoved)
time=np.delete(time,toBeRemoved)

thickness_err = 0.5*np.ones(thickness.shape)

count_err = np.empty(count.shape)
#switch over to cm thicknesses
thickness = 0.1*thickness

for (i,c) in enumerate(count):
    count_err[i] = np.sqrt(c)

time_err = np.ones(time.shape)

intensity = count/time
# remove sqrt if errors are too small
intensity_err = intensity*np.sqrt((count_err/count)**2+(time_err/time)**2)

#partition in terms of energy
energydata = dict()
for i in set(energy):
    energydata[i]=([],[],[])

for i,t in enumerate(thickness):
    energydata[energy[i]][0].append(t)

for i,t in enumerate(intensity):
    energydata[energy[i]][1].append(t)

for i,t in enumerate(intensity_err):
    energydata[energy[i]][2].append(t)



class Parameter:
    def __init__(self, initialvalue, name):
        self.value = initialvalue
        self.name=name
    def set(self, value):
        self.value = value
    def __call__(self):
        return self.value

def fit(function, parameters, x, data, u):
    def fitfun(params):
        for i,p in enumerate(parameters):
            p.set(params[i])
        return (data - function(x))/u
 
    if x is None: x = arange(data.shape[0])
    if u is None: u = ones(data.shape[0],"float")
    p = [param() for param in parameters]
    return leastsq(fitfun, p, full_output=1)

info = dict()
att_coef = dict()
att_errs = dict()

iron_const = 8.491767433789953e+22*26
for e in energydata.keys():
    A = Parameter(1.,'A')
    B = Parameter(1.,'B')
    C = Parameter(1.,'C')
    p0 = [A,B,C]
    def expfit(r):
        return A()*exp(-r*B())+C()
    (p2,cov,info,msg,success) = fit(expfit,p0,np.array(energydata[e][0]), np.array(energydata[e][1]), np.array(energydata[e][2]))
    chisq = sum(info['fvec']*info['fvec'])
    dof = len(energydata[e][0])-len(p0)
    print "Energy ",e,A(),B(),B()/iron_const
    att_coef[e] = B()
    att_errs[e] = sqrt(cov[1,1])*sqrt(chisq/dof)
    print "Fitted parameters at minimum, with 68% C.I.:"
    for i,pmin in enumerate(p2):
        print "%2i %-10s %12f +/- %10f"%(i,p0[i].name,pmin,sqrt(cov[i,i])*sqrt(chisq/dof))
    print

    fig = matplotlib.pyplot.figure()
    errorbar(np.array(energydata[e][0]),np.array(energydata[e][1]),yerr=np.array(energydata[e][2]),xerr=0.005,fmt='r.',label='Data')
    xplot = linspace(-0.1,6,70)
    yplot = expfit(xplot)
    ax = axes()
    ax.plot(xplot,yplot,'b-',label='Fit')
    ax.text(0.4,0.7,"Fit equation: $%.3fe^{-%.3fx}+(%.3f)$\nReduced Chi-Sq: %.3f\n"%(A(),B(),C(),chisq/dof),fontsize=12,
        horizontalalignment='left',transform = ax.transAxes)
    xlabel('Absorber thickness(cm)')
    ylabel('Intensity(counts/second)')
    title(str(e)+'keV data and fit')
    legend()
    fig.savefig(str(e)+'keV.png',format='png')
    
figure()
ax = axes()
es = np.array(energydata.keys())
cs = np.empty(es.shape)
ce = np.empty(es.shape)
for (i,e) in enumerate(es):
    cs[i] = att_coef[e]
    ce[i] = att_errs[e]

ax.errorbar(es,cs,yerr=ce,fmt='r.',label='Attenuation coefficients')
ax.set_xscale('log')
ax.set_yscale('log')
xlabel('log(energy)')
ylabel('log(attenuation coefficient)')
savefig('att_coefs.png',format='png')
for i in range(len(energydata.keys())):
    print '$'+str(es[i])+'$ & $'+str(cs[i])+'\\pm'+str(ce[i])+'$ & $' + str(cs[i]/iron_const)+"\\pm"+str(ce[i]/iron_const)+'$\\\\'

#now fit the gaussian to the peak
chan,chanct=loadtxt('1330peak',unpack=True,usecols=[0,1])
cterr = sqrt(chanct)+1e-5

figure()
title('The 1330keV peak of Na-22 with a Gaussian fit')
R = Parameter(3500.,'R')
u = Parameter(480.,'u')
s = Parameter(6,'s')
p0=[R,u,s]
def gaussian(x):
    return R()*exp(-((x-u())**2)/(s()**2))
p2,cov,info,mesg,success=fit(gaussian, p0, chan, chanct, cterr)

chisq = sum(info['fvec']*info['fvec'])
dof = len(chan)-len(p0)

freqs = linspace(445,505,200)
fittedplot = gaussian(freqs)

errorbar(chan,chanct,yerr=cterr,fmt='r.',label='Count data')
ax=axes()
ax.plot(freqs,fittedplot,'b-',label='fit')
ax.text(0.4,0.2,"Mean: $%.3f\\pm%.3f$\nStandard deviation:$%.3f$"%(u(),sqrt(cov[1,1]*sqrt(chisq/dof)),s()),fontsize=12,
    horizontalalignment='left',transform = ax.transAxes)
xlabel('Channel')
ylabel('Counts')
legend()
savefig('1330peakplot.png',format='png')
#show()

