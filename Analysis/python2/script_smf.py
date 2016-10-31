from matplotlib.pyplot import *

# Haven't worked out why but show() does not replot the data.
# So need to close the plotting window between calls.

# Parameters for Hen14
hubble_Hen14=0.673
boxside_Hen14=480.28  # Units Mpc/h 
datafile_Hen14='data/Hen14_sfh2/gal_Hen14.pkl'

# Parameters for Guo10
hubble_Guo10=0.73
boxside_Guo10=500.  # Units Mpc/h 
datafile_Guo10='data/Guo10_sfh2/gal_Guo10.pkl'

# Define limits of plot
# Because I want to use these as plotting data, they need to be arrays, 
# not lists.  The following seems to work
binperdex=10
xrange=np.array([7,12])
nbin=(xrange[1]-xrange[0])*binperdex

#--------------------------------------------------------------------

# Loading in the galaxy data is slow, so we only want to do it once.
# This checks to see whether or not we have done so.
# To force a reread do "del gals"
# Load in Hen14
try:
    mstar_Hen14
except:
    pickleFile=datafile_Hen14
    execfile('script_unpickle.py')
    # Mass in units of Msun/h^2
    mstar_Hen14=(gals['DiskMass']+gals['BulgeMass'])*1e10*hubble_Hen14
# Load in Guo10
try:
    mstar_Guo10
except:
    pickleFile=datafile_Guo10
    execfile('script_unpickle.py')
    # Mass in units of Msun/h^2, converting to the Planck h
    mstar_Guo10=(gals['DiskMass']+gals['BulgeMass'])*1e10 \
        *hubble_Hen14**2/hubble_Guo10

# Put into bins and normalise to number per unit volume (Mpc/h) per dex
#Hen14
nobj,bins,junk=hist(np.log10(mstar_Hen14), bins=nbin, range=xrange, log=True)
y_Hen14=nobj/(boxside_Hen14)**3*binperdex
#Guo10
nobj,bins,junk=hist(np.log10(mstar_Guo10), bins=nbin, range=xrange, log=True)
y_Guo10=nobj/(boxside_Guo10)**3*binperdex

# Plot at centre of bins
x=0.5*(bins[:-1]+bins[1:])

# Plot
close()
semilogy(x,y_Hen14,'+r',x,y_Guo10,'*b')
axis([12.4,7,10.**(-5.9),10.**0.5])
xlabel(r'$\log_{10}(M_*/h^{-2}M_\odot)$')
ylabel(r'$\log_{10}(N/(\mathrm{dex}\ (h^{-1}\mathrm{Mpc})^3)$')
legend(['Hen14','Guo10'],2)
grid(True)
show()
savefig('smf.png')

