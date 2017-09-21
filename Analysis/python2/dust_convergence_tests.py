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



#---------------------All NEW dust Plots
#fig = plt.figure(figsize=(9,9))
fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(9,9))
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)

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
	
	if(loop==1 or loop==4 or loop==2 or loop==3):
		plt.tick_params(axis='both', which='major', labelsize=12,width=2,length=6, labelleft ='off', labelbottom='off', direction='inout')
		plt.tick_params(axis='both', which='minor', labelsize=12,width=2,length=6, labelleft ='off', labelbottom='off', direction='inout')
	
	if(loop==5 or loop==6 or loop==7):
		plt.tick_params(axis='both', which='major', labelsize=12,width=2,length=6, labelleft ='off', labelbottom='on', direction='inout')
		plt.tick_params(axis='both', which='minor', labelsize=12,width=2,length=6, labelleft ='off', labelbottom='on', direction='inout')
		plt.xlim([6.0,11.98])
		plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
	
	if(loop==8):
		plt.tick_params(axis='both', which='major', labelsize=12,width=2,length=6, labelleft ='off', labelbottom='on', direction='inout')
		plt.tick_params(axis='both', which='minor', labelsize=12,width=2,length=6, labelleft ='off', labelbottom='on', direction='inout')
		plt.xlim([6.01,12.00])
		plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)

	if(loop==0 or loop==3):
		plt.tick_params(axis='both', which='major', labelsize=12,width=2,length=6, labelleft ='on', labelbottom='off', direction='inout')
		plt.tick_params(axis='both', which='minor', labelsize=12,width=2,length=6, labelleft ='on', labelbottom='off', direction='inout')
		plt.ylim([0.01,10])
		plt.ylabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)

	if(loop==6):
		plt.tick_params(axis='both', which='major', labelsize=12,width=2,length=6, labelleft ='on', labelbottom='on', direction='inout')
		plt.tick_params(axis='both', which='minor', labelsize=12,width=2,length=6, labelleft ='on', labelbottom='on', direction='inout')
		plt.ylim([0.01,10])
		plt.xlim([6.01,11.98])
		plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
		plt.ylabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)

	plt.tick_params(axis='both', which='major', labelsize=12,width=2,length=6, direction='inout')
	plt.tick_params(axis='both', which='minor', labelsize=12,width=2,length=6, direction='inout')
	
	plt.plot([9.0,9.0],[-1.0,11.0],linewidth=2,color='b')

axes = fig.get_axes()

for ax in axes:
	[i.set_linewidth(2.1) for i in ax.spines.itervalues()]

pylab.savefig('./graphs/convergence_test_ALL_.png', bbox_inches=0)
plt.close()
	










