#############################################################################
#   Fit Absorption Curve Data
#
#   This script reads a user generated text file containing Gamma intensities
#   vs. absorber thickness data and then performs a fit on the data to determine
#   the interaction cross section.
#
#   The fitting routine used is the Levenberg-Marquardt method. This script is
#   more complicated then usual, but most of the complexity can be ignored by
#   student.
#
#   Modified 101912 to print out the uncertainties of the fit parameters.
#       MCC
#############################################################################

# Import required libraries
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
 
def fit(function, parameters, x, data, u):
    def fitfun(params):
        for i,p in enumerate(parameters):
            p.set(params[i])
        return (data - function(x))/u
 
    if x is None: x = arange(data.shape[0])
    if u is None: u = ones(data.shape[0],"float")
    p = [param() for param in parameters]
    return leastsq(fitfun, p, full_output=1)
#
#############################################################################

# Read in the data to be fit from a user created text file.  The script
# assumes three columns of data.
# Column 1 = Absorber thicknesses.
# Column 2 = Gamma intensities.
# Column 3 = Uncertainties for Gamma intensities.
x, y, ye =loadtxt('abs_data.dat',unpack=True, usecols=[0,1,2])

# Generate figure object.
fig = matplotlib.pyplot.figure()

# Plot the data with error bars.
errorbar(x,y, yerr=ye, fmt='', label="")

# Define a simple exponential decay function to be fitted.
# The free parameters of the fit are:
# A = The amplitude of the curve.
# B = The time constant of the decay.
def Expdecay(fr):
    return A()*exp(-fr*B())
            
# Define the fit patameters and assign them initial guess values.
A=Parameter(1.0,"A")
B=Parameter(1.0,"B")

# Put the fit patameters into a list.
p0=[A,B]

# Do fit using Levenberg-Marquardt method.  We provide the fitting routine
# with:
# Expdecay = The function to be fit as defined above.
# p0 = The list of fit parameters and their initial values as defined above.
# x = The x data values as read in from the user data file.
# y = The y data values as read in from the user data file.
# y = The uncertainties in the y data values as read in from the user data file.
#
# The fitting routine returns information about the success of the fit and the
# final best fit values for the parameters.  The relevent information returned
# is:
# success = An integer indicating whether or not the fit successfully
#           converged on a solution.  A value of 1 indicates success.
# mesg = A descriptive message of why the fit succeded or failed.
# info = Information about the goodness of the fit.
# p2 and cov = Other information which we do not make use of.
#
p2,cov,info,mesg,success=fit(Expdecay, p0, x, y, ye)

# If the fit succeeded, as indicated by the value of 'success' being 1, print
# the message "Converged" to the terminal.  If the fit did not succeed print
# "Not converged" along with the contents of the 'mesg' returned by the fitting
# routine.
if success==1:
    print "Converged"
else:
    print "Not converged"
    print mesg

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

# Generate 50 x axis values over the range of absorber thicknesses used in
# the experimental data.
thick=linspace(0.0,5,50)

# Using the exponential decay function defined earlier, along with the best
# fit values for th free parameters and the x values generated above, calculate
# y values.
initialplot=Expdecay(thick)

# Create a new axes object within which the best fit data will be plotted.
ax = axes()

# Plot the best fit data.
ax.plot(thick, Expdecay(thick), '', label="")

# Display the best fit parameters on the plot.
# The location of the text may need to be adjusted.
ax.text(0.4, 0.6,
        '$Best$ $Fit$ $Values$ : $A*exp(-x*B)$ \n $A$ =  %.3f \n $B$ =  %.3f \n $Reduced$ $CHISQ$ =  %.3f ' % (A(),B(),chisq/dof),
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

