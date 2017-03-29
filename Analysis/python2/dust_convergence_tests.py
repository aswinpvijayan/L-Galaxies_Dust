# Script to read in pickled file of galaxy properties.
# Requires pickleFile to hold name of desired input file.
# Returns galaxy properties in gals.

def fit_scatter(x, y, ret_n=False, ret_sterr=False, nbins=10):
    '''
    Bins scattered points and fits with error bars
    '''
    import numpy as np

    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)

    bin_centres = (_[1:] + _[:-1])/2.

    if ret_sterr:
        stderr = std/np.sqrt(n)
        if ret_n:
            return bin_centres, mean, std, stderr, n
        return bin_centres, mean, std, stderr

    if ret_n:
        return bin_centres, mean, std, n
    return bin_centres, mean, std
    
    
    


import cPickle
import numpy as np
import matplotlib.pyplot as plt
import pylab
import sys


#### Read in observational data

#### Remy-Ruyer 2014 z=0
#z0_obs_DM, obs_DM_err, obs_SM = np.loadtxt('../observations/z0.txt',unpack=True)
#z0_obs_DM_act_err = obs_DM*(obs_DM_err/100.0)
#Name DustMass DMerror% SMass MHI MHIerror% 12+log(O/H) MH2,mw MH2,Z
Remy_DM, Remy_DMerr, Remy_SM, Remy_MHI, Remy_MHIerr, Remy_Metals, Remy_H2mw, Remy_MH2z = np.loadtxt('../observations/Remy_Ruyer_2014_KINGFISH_z0.txt',unpack=True,comments='#')
Remy_DM_err_actual = Remy_DM*(Remy_DMerr/100.0)
Remy_Dust_Gas_Ratio = Remy_DM / (Remy_MHI + Remy_H2mw)
Remy_Dust_Metal_Ratio = Remy_DM / (((10**(Remy_Metals - 8.69))*0.02)*(Remy_MHI + Remy_H2mw))

#print Remy_SM, np.log10(Remy_Dust_Metal_Ratio), Remy_DM, (((10**(Remy_Metals - 8.69))*0.02)*(Remy_MHI + Remy_H2mw))


#### Santini 2014 z=0,1,2
# SM DM DM+err DM-err
Santini_SM_z0, Santini_DM_z0, Santini_DMuperr_z0, Santini_DMdownerr_z0 = np.loadtxt('../observations/Santini_2014_z0.txt',unpack=True,comments='#')
Santini_SM_z1, Santini_DM_z1, Santini_DMuperr_z1, Santini_DMdownerr_z1 = np.loadtxt('../observations/Santini_2014_z1.txt',unpack=True,comments='#')
Santini_SM_z2, Santini_DM_z2, Santini_DMuperr_z2, Santini_DMdownerr_z2 = np.loadtxt('../observations/Santini_2014_z2.txt',unpack=True,comments='#')

#### daCunha 2015
#Phot_z z+err z-err SM SM+err SM-err DM DM+err DM-err
daCunha_z,daCunha_zuperr,daCunha_zdownerr,daCunha_SM,daCunha_SMuperr,daCunha_SMdownerr,daCunha_DM,daCunha_DMuperr,daCunha_DMdownerr = np.loadtxt('../observations/daCunha_2015_z_3_4_5.txt',unpack=True,comments='#')

daCunha_z2,daCunha_zuperr_z2,daCunha_zdownerr_z2,daCunha_SM_z2,daCunha_SMuperr_z2,daCunha_SMdownerr_z2,daCunha_DM_z2,daCunha_DMuperr_z2,daCunha_DMdownerr_z2 = np.loadtxt('../observations/daCunha_2015_z_2.txt',unpack=True,comments='#')
daCunha_z3,daCunha_zuperr_z3,daCunha_zdownerr_z3,daCunha_SM_z3,daCunha_SMuperr_z3,daCunha_SMdownerr_z3,daCunha_DM_z3,daCunha_DMuperr_z3,daCunha_DMdownerr_z3 = np.loadtxt('../observations/daCunha_2015_z_3.txt',unpack=True,comments='#')
daCunha_z4,daCunha_zuperr_z4,daCunha_zdownerr_z4,daCunha_SM_z4,daCunha_SMuperr_z4,daCunha_SMdownerr_z4,daCunha_DM_z4,daCunha_DMuperr_z4,daCunha_DMdownerr_z4 = np.loadtxt('../observations/daCunha_2015_z_4.txt',unpack=True,comments='#')
daCunha_z5,daCunha_zuperr_z5,daCunha_zdownerr_z5,daCunha_SM_z5,daCunha_SMuperr_z5,daCunha_SMdownerr_z5,daCunha_DM_z5,daCunha_DMuperr_z5,daCunha_DMdownerr_z5 = np.loadtxt('../observations/daCunha_2015_z_5.txt',unpack=True,comments='#')
daCunha_z6,daCunha_zuperr_z6,daCunha_zdownerr_z6,daCunha_SM_z6,daCunha_SMuperr_z6,daCunha_SMdownerr_z6,daCunha_DM_z6,daCunha_DMuperr_z6,daCunha_DMdownerr_z6 = np.loadtxt('../observations/daCunha_2015_z_6.txt',unpack=True,comments='#')


####Mancini2015

Mancini_z, Mancini_SM, Mancini_SMerr, Mancini_DM, Mancini_DMerr = np.loadtxt('../observations/Mancini_2015_z6_z7.txt',unpack=True,comments='#')

###Bourne2012

Bourne_MEDZ,Bourne_MEDM,Bourne_MEDC,Bourne_MLIMITS_down,Bourne_MLIMITS_up,Bourne_MEDMSTAR,Bourne_MEDMSTARERR,Bourne_MEDMDUST,Bourne_MEDMDUSTERR,Bourne_MEDDMS,Bourne_MEDDMSERR,Bourne_NBIN = np.loadtxt('../observations/Bourne2012_z0.txt',unpack=True,comments='#')


# print 'Number of arguments:', len(sys.argv), 'arguments.'
# print 'Argument List:', str(sys.argv)

#---------------------All NEW dust Plots
fig = plt.figure(figsize=(9,9))
for loop in range(0,9):

	SM_bins_MR,Dust_bins_MR,Dust_std_dev_MR,Dust_std_err_MR,count_MR = np.loadtxt('./binned_data/SM_DM_z'+str(loop)+'_MR.txt',unpack=True)
	SM_bins_MRII,Dust_bins_MRII,Dust_std_dev_MRII,Dust_std_err_MRII,count_MRII = np.loadtxt('./binned_data/SM_DM_z'+str(loop)+'_MRII.txt',unpack=True)
	
	plt.subplot(3,3,loop+1)
	
	plt.xlim([6,12])
	plt.ylim([0,10])
	
	plt.errorbar(SM_bins_MR,Dust_bins_MR,yerr=(Dust_std_err_MR),color='k',label='MR')
	plt.errorbar(SM_bins_MRII,Dust_bins_MRII,yerr=(Dust_std_err_MRII),color='r',label='MRII')
	if loop==7:
		plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
	if loop==3:
		plt.ylabel(r'log$_{10}$(M$_{\rm{d}}$/M$_{\odot}$)', fontsize=18)
	
	plt.text(10,1,"z="+str(loop),fontsize=18)
	
	plt.tick_params(axis='both', which='major', labelsize=10,width=2)
	plt.tick_params(axis='both', which='minor', labelsize=10,width=2)

axes = fig.get_axes()

for ax in axes:
	[i.set_linewidth(2.1) for i in ax.spines.itervalues()]

pylab.savefig('./graphs/convergence_test_ALL_.png', bbox_inches=0)
plt.close()
	










