#!/usr/bin/python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from sig_fig import *

mag,a1,da1,a2,da2,b1,db1,b2,db2,c1,dc1,c2,dc2=loadtxt('peakdat', unpack=True)
#calculate actual radii and errors
a_rad = (a1-a2)*0.5
b_rad = (b1-b2)*0.5
c_rad = (c1-c2)*0.5
da_rad = np.sqrt(da1**2+da2**2)/2
db_rad = np.sqrt(db1**2+db2**2)/2
dc_rad = np.sqrt(dc1**2+dc2**2)/2

#put stuff in dicts
data = dict()
for b in set(mag):
    data[b]=[[],[],[],[],[],[]]

for i in range(len(mag)):
    data[mag[i]][0].append(a_rad[i])
    data[mag[i]][1].append(da_rad[i])
    data[mag[i]][2].append(b_rad[i])
    data[mag[i]][3].append(db_rad[i])
    data[mag[i]][4].append(c_rad[i])
    data[mag[i]][5].append(dc_rad[i])
#import the .103V data
ra,rb,rc,dra,drb,drc=loadtxt('peakdat103',unpack=True)
data[103.]=[ra,dra,rb,drb,rc,drc]

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

epsilons = dict()
for b in data.keys():
    epsilons[b]=[0,0,0,0,0,0]

ringdict = {0:'a',2:'b',4:'c'}
fpspace = 0.99
dfpspace = 0.001
for b in data.keys():
    figure()
    axes()
    formatkey = {0:('r.','g-'),2:('b.','y-'),4:('c.','m-')}
    add_txt=''
    for r in [0,2,4]:
#        print len(data[b][r+1])
#        print np.arange(1,len(data[b][r])+1)
#        print np.array(data[b][r])**2
#        print np.array(data[b][r+1])*2
        A = Parameter(1.,'A')
        B = Parameter(1.,'B')
        if b==225 and r==4:
            B.set(10000.)
        p0 = [A,B]
        def lin_fit(x):
            return A()*x+B()
      #  if b==225 and r==4:
       #     for i in range(len(data[b][r+1])): data[b][r+1][i]=0.01
        p2,cov,info,msg,success = fit(lin_fit,p0,np.arange(1,len(data[b][r])+1),np.array(data[b][r])**2,np.array(data[b][r+1])*np.sqrt(2)*np.array(data[b][r]))
        chisq = sum(info['fvec']*info['fvec'])
        dof = len(data[b][r])-len(p0)
        epsilons[b][r]=B()/A()+1
        a_err = np.sqrt(cov[0,0])*np.sqrt(chisq/dof)
        b_err = np.sqrt(cov[1,1])*np.sqrt(chisq/dof)
        frac_err = np.sqrt((b_err/B())**2+(a_err/A())**2)
        epsilons[b][r+1]=frac_err*epsilons[b][r]
        print "Magnetic field: %d, ring %d"%(b,r)
        print "Converged with chi squared ",chisq
        print "degrees of freedom, dof ", dof
        print "RMS of residuals (i.e. sqrt(chisq/dof)) ", sqrt(chisq/dof)
        print "Reduced chisq (i.e. variance of residuals) ", chisq/dof
        print "a, da: ",A(), a_err
        print "b, db: ",B(), b_err
        print "epsilon: "+sig_fig(epsilons[b][r],epsilons[b][r+1])
        print 
        
#B-field calibration
volt, bmtr, dbmtr = loadtxt('calibrationvb',unpack=True)
bref,bmtr1,dbmtr1 = loadtxt('calibrationrm',unpack=True)

A=Parameter(1.,'A')
B=Parameter(1.,'B')
a=Parameter(1.,'a')
b=Parameter(1.,'b')
p0=[A,B]
p1=[a,b]
p00,cov0,info0,msg0,success0=fit((lambda x:A()*x+B()),p0,volt,bmtr,dbmtr)
p11,cov1,info1,msg1,success1=fit((lambda x:a()*x+b()),p1,bref,bmtr1,dbmtr1)
dof0=len(volt)-2
dof1=len(bref)-2
chisq0=sum(info0['fvec']**2)
chisq1=sum(info1['fvec']**2)
print "First fit reduced chisq: %.3f, second fit: %.3f"%((chisq0/dof0),(chisq1/dof1))
A_err=np.sqrt(cov0[0,0])*np.sqrt(chisq0/dof0)
B_err=np.sqrt(cov0[1,1])*np.sqrt(chisq0/dof0)
a_err=np.sqrt(cov1[0,0])*np.sqrt(chisq1/dof1)
b_err=np.sqrt(cov1[1,1])*np.sqrt(chisq1/dof1)

slope=A()/a()
dslope=slope*np.sqrt((A_err/A())**2+(a_err/a())**2)

yint=(B()-b())/a()
dyint=np.sqrt((B_err/a())**2+(b_err/a())**2+((B()-b())*a_err/a()**2)**2)
#print "B_actual=%sv+%s"%(sig_fig(slope,dslope),sig_fig(yint,dyint))
mag1 = np.concatenate((np.array([103]),sort(list(set(mag)))))
b_actuals = slope*mag1/1000+yint
db_actuals = np.sqrt((slope*mag1/1000)**2*((dslope/slope)**2+(1/mag1)**2)+dyint**2)

bconvs = dict()
for i in range(len(mag1)):
    bconvs[mag1[i]]=(b_actuals[i],db_actuals[i])

print b_actuals
print db_actuals
fudge=1.23984e-4*1000000
magneton=0.002894191

bb = []
dbb = []
eb = []
deb = []
for b in sort(epsilons.keys()):
    print "Voltage: %.3f corresponds to %.3f +/- %.3fT"%(b,bconvs[b][0],bconvs[b][1])
    print "a-b difference: "+sig_fig((epsilons[b][0]-epsilons[b][2])/(2*fpspace),np.sqrt((0.5*epsilons[b][1]/fpspace)**2+(0.5*epsilons[b][3]/fpspace)**2+((epsilons[b][0]-epsilons[b][2])/(2*fpspace**2)*dfpspace)**2))
    print "a-c difference: "+sig_fig((epsilons[b][0]-epsilons[b][4])/(2*fpspace),np.sqrt((0.5*epsilons[b][1]/fpspace)**2+(0.5*epsilons[b][5]/fpspace)**2+((epsilons[b][0]-epsilons[b][4])/(2*fpspace**2)*dfpspace)**2))
    print "b-c difference: "+sig_fig((epsilons[b][2]-epsilons[b][4])/(2*fpspace),np.sqrt((0.5*epsilons[b][3]/fpspace)**2+(0.5*epsilons[b][5]/fpspace)**2+((epsilons[b][2]-epsilons[b][4])/(2*fpspace**2)*dfpspace)**2))
    text = "%.3f & $%s$"%(b/1000,sig_fig(bconvs[b][0],bconvs[b][1]))
    text = text + " & $" + sig_fig(-fudge*(epsilons[b][0]-epsilons[b][2])/(2*fpspace),fudge*np.sqrt((0.5*epsilons[b][1]/fpspace)**2+(0.5*epsilons[b][3]/fpspace)**2+((epsilons[b][0]-epsilons[b][2])/(2*fpspace**2)*dfpspace)**2))
    text = text + "$ & $" + sig_fig(-fudge*(epsilons[b][2]-epsilons[b][4])/(2*fpspace),fudge*np.sqrt((0.5*epsilons[b][3]/fpspace)**2+(0.5*epsilons[b][5]/fpspace)**2+((epsilons[b][2]-epsilons[b][4])/(2*fpspace**2)*dfpspace)**2))
    text = text + "$ & $" + sig_fig(magneton*bconvs[b][0],magneton*bconvs[b][1])
    text = text + "$\\\\"
    print text
    bb.append(bconvs[b][0])
    dbb.append(bconvs[b][1])
    eb.append(-fudge*(epsilons[b][2]-epsilons[b][4])/(2*fpspace))
    deb.append(fudge*np.sqrt((0.5*epsilons[b][3]/fpspace)**2+(0.5*epsilons[b][5]/fpspace)**2+((epsilons[b][2]-epsilons[b][4])/(2*fpspace**2)*dfpspace)**2))


S = Parameter(1.,'S')
T = Parameter(1.,'T')
p0 = [S,T]
def lin_fit(x):
    return S()*x+T()
p2,cov,info,msg,success = fit(lin_fit,p0,np.array(bb),np.array(eb),np.array(deb))
chisq = sum(info['fvec']**2)
print 'Reduced chisq of magneton fit: %.3f'%(chisq/6)
print 'Bohr magneton: %s'%sig_fig(S()*2,2*np.sqrt(cov[0,0])*np.sqrt(chisq/6))
print '('+sig_fig(S(),np.sqrt(cov[0,0])*np.sqrt(chisq/6))+')\\cdot B + ('+sig_fig(T(),np.sqrt(cov[1,1])*np.sqrt(chisq/6))+')'

figure()
errorbar(np.array(bb),np.array(eb),fmt='k.',yerr=np.array(deb),xerr=np.array(dbb))
indepvar = range(1,9)
rg = np.arange(0,10000,10)
plot(rg,lin_fit(rg),'b-')
gca().set_xlabel('B-field (gauss)')
gca().set_ylabel('Energy splitting (micro electron volts)')
savefig('bc.png')
#how()
#@        if chisq/dof>1:
#@            print b,r
#@            errorbar(np.arange(1,len(data[b][r])+1),np.array(data[b][r])**2,fmt='b+',yerr=np.array(data[b][r+1])*2)
#@            xs = arange(10)
#@            plot(xs,lin_fit(xs),'k-')
#@            show()

tables = open('tables','w')
for b in data.keys():
    tables.write('\\begin{tabular}{|c|c|}\n')
    tables.write('\\hline\n')
    tables.write('\\multicolumn{2}{|c|}{Radii for V=%.3f}\\\\\n'%b)
    tables.write('\\hline\n')
    for r in [0,2,4]:
        tables.write('Ring '+ringdict[r]+' & $')
        for i in range(len(data[b][r])):
            tables.write(sig_fig(data[b][r][i],data[b][r+1][i])+', ')
        tables.write('$\\\\\n')
    tables.write('\\hline\n')
    tables.write('\\end{tabular}')
    tables.write('\n\n')

plottex = open('plottex','w')
for b in data.keys():
    plottex.write('\\begin{figure}\n')
    plottex.write('\\centering\n')
    plottex.write('\\includegraphics[width=500pt]{dm0/plots/%d-rplot.png}\n'%b)
    plottex.write('\\end{figure}\n')
    plottex.write('\\vspace{10pt}\n\n')
