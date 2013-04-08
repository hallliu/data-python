# Python example for analysis of Electron Spin Resonance experiment
# Includes uncertainties in x in the chi squared fit
# P22101, Autumn 2012
# 11/19/2012

# import required libraries
from pylab import *
from scipy.optimize import leastsq
from matplotlib.backends.backend_pdf import PdfPages
import numpy as numpy

#############################################################################
#   From  www.scipy.org/Cookbook/FittingData :
#   Modified to add uncertainties, full output
#
#   The following section contains definitions necessary for the curve
#   fit method imployed.  Students with programming experience are
#   encouraged to read through the definitions to understand what they
#   are doing.  However this section can be safely ignored by students
#   unfamiliar with the concepts of classes and parameters.
#
class Parameter:
    def __init__(self, initialvalue, name):
        self.value = initialvalue
        self.name=name
    def set(self, value):
        self.value = value
    def __call__(self):
        return self.value
 
def fit(function, parameters, x, data, v, u):
    def fitfun(params):
        for i,p in enumerate(parameters):
            p.set(params[i])
        return (data - function(x))/sqrt(u**2+(v*A())**2)
# modified to also take errors in x into account. Note the function now takes
# an additional argument, v, the error in x
 
    if x is None: x = arange(data.shape[0])
    if u is None: u = ones(data.shape[0],"float")
    p = [param() for param in parameters]
    return leastsq(fitfun, p, full_output=1)
# ========================================

# The following lines have example (and mostly non-sensical) data. You
# should substitute these for your own results.

x=array([1,2,3,4,5])
xe=array([.1,.1,.1,.1,.1])
y=array([.9,1.9,3.3,4.5,5.2])
ye=array([.1,.1,.1,.1,.1])

# End of data input.

# generate figure object
fig = matplotlib.pyplot.figure()

# plot the data with error bars
errorbar(x,y, yerr=ye, xerr=xe, fmt='', label="")

# define function to be fitted: a simple line, with parameters m (slope) and
# b (offset) from 0 at x=0
def line(x):
    return A()*x+B()
            

# Define the fit parameters and assign them initial guess values.
# Note that with uncertainty in x, it becomes much more important to get this
# guess close to the real values or the fit may not work. You may need to try
# more than once.
A=Parameter(1.0,"A")
B=Parameter(1.0,"B")

# Put the fit parameters into a list.
p0=[A,B]

# Do fit using Levenberg-Marquardt method.  We provide the fitting routine
# with the usual stuff, plus xe, the uncertainty in x. Note the ordering
# of the arguments is important.
p2,cov,info,mesg,success=fit(line, p0, x, y, xe, ye)

# Using the data returned by the fitting routine in "info", calculate the
# final chi square value for the fit.
chisq=sum(info["fvec"]*info["fvec"])

# Determine the number of degrees of freedom for the fit.
dof=len(x)-len(p0)

# Print data about the goodness of the fit to the terminal.
print "Converged with chi squared ",chisq
print "degrees of freedom, dof ", dof
print "RMS of residuals (i.e. sqrt(chisq/dof)) ", sqrt(chisq/dof)
print "Reduced chisq (i.e. variance of residuals) ", chisq/dof
print

# uncertainties are calculated as per gnuplot, "fixing" the result
# for non unit values of the reduced chisq.
# values at min match gnuplot
print "Fitted parameters at minimum, with 68% C.I.:"
for i,pmin in enumerate(p2):
    print "%2i %-10s %12f +/- %10f"%(i,p0[i].name,pmin,sqrt(cov[i,i])*sqrt(chisq/dof))
print

print "Correlation matrix"
# correlation matrix close to gnuplot
print "               ",
for i in range(len(p0)): print "%-10s"%(p0[i].name,),
print
for i in range(len(p2)):
    print "%10s"%p0[i].name,
    for j in range(i+1):
        print "%10f"%(cov[i,j]/sqrt(cov[i,i]*cov[j,j]),),
    print

#############################################################################
# Now we wish to plot the best fit on the same graph as the plot of the data.
# To do so we will use the best fit parameters returned by the fitting routine
# to calculate y values corresponding to the same range of x values covered
# by the data.
#

# Generate 50 x axis values over the range of x used in
# the experimental data.
pts=linspace(0.0,x.max(),50)

# Using the fit function defined earlier, along with the best
# fit values for the free parameters and the x values generated above, calculate
# y values.
initialplot=line(pts)

# Create a new axes object within which the best fit data will be plotted.
ax = axes()

# Plot the best fit data.
ax.plot(pts, line(pts), '', label="")

# Display the best fit parameters on the plot.
# The location of the text may need to be adjusted.
ax.text(0.4, 0.6,
        '$Best$ $Fit$ $Values$ : $A*x+B$ \n $A$ =  %.3f \n $B$ =  %.3f \n $Reduced$ $CHISQ$ =  %.3f ' % (A(),B(),chisq/dof),
        fontsize=12,
        horizontalalignment='left',
        verticalalignment='center',
        transform=ax.transAxes)

# Add your own axis labels and legend.
xlabel("")
ylabel("")

legend("")

fig.set_facecolor('white')

# Save the plot to a PDF file.
plt.savefig('.pdf', bbox_inches=0)

# Show the durn plot already!
show()
 
