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




print "Redshift [Number of galaxies in each mass bin]"	

#for loop in range(0,10):
for loop in range(0,1):

	fin = open('../data/metals/lgal_z'+str(loop)+'.pkl','rb')
	gals_metals=cPickle.load(fin)
	fin.close()

	fin = open('../data/lgal_z'+str(loop)+'.pkl','rb')
	gals_dust=cPickle.load(fin)
	fin.close()

	Stellar_Mass_metals = np.zeros(len(gals_metals['Type']))
	SFR_metals = np.zeros(len(gals_metals['Type']))
	
	H_Mass_metals = np.zeros(len(gals_metals['Type']))
	He_Mass_metals = np.zeros(len(gals_metals['Type']))
	Cb_Mass_metals = np.zeros(len(gals_metals['Type']))
	N_Mass_metals = np.zeros(len(gals_metals['Type']))
	O_Mass_metals = np.zeros(len(gals_metals['Type']))
	Ne_Mass_metals = np.zeros(len(gals_metals['Type']))
	Mg_Mass_metals = np.zeros(len(gals_metals['Type']))
	Si_Mass_metals = np.zeros(len(gals_metals['Type']))
	S_Mass_metals = np.zeros(len(gals_metals['Type']))
	Ca_Mass_metals = np.zeros(len(gals_metals['Type']))
	Fe_Mass_metals = np.zeros(len(gals_metals['Type']))
	
	for i in range(0,len(gals_metals['Type'])):
		H_Mass_metals[i] = gals_metals['ColdGas_elements'][i][0]
		He_Mass_metals[i] = gals_metals['ColdGas_elements'][i][1]
		Cb_Mass_metals[i] = gals_metals['ColdGas_elements'][i][2]
		N_Mass_metals[i] = gals_metals['ColdGas_elements'][i][3]
		O_Mass_metals[i] = gals_metals['ColdGas_elements'][i][4]
		Ne_Mass_metals[i] = gals_metals['ColdGas_elements'][i][5]
		Mg_Mass_metals[i] = gals_metals['ColdGas_elements'][i][6]
		Si_Mass_metals[i] = gals_metals['ColdGas_elements'][i][7]
		S_Mass_metals[i] = gals_metals['ColdGas_elements'][i][8]
		Ca_Mass_metals[i] = gals_metals['ColdGas_elements'][i][9]
		Fe_Mass_metals[i] = gals_metals['ColdGas_elements'][i][10]
		
		Stellar_Mass_metals[i] = gals_metals['StellarMass'][i]*1.0E10/0.673
		SFR_metals[i] = gals_metals['Sfr'][i]

	Stellar_Mass_dust = np.zeros(len(gals_dust['Type']))
	SFR_dust = np.zeros(len(gals_dust['Type']))

	H_Mass_dust = np.zeros(len(gals_dust['Type']))
	He_Mass_dust = np.zeros(len(gals_dust['Type']))
	Cb_Mass_dust = np.zeros(len(gals_dust['Type']))
	N_Mass_dust = np.zeros(len(gals_dust['Type']))
	O_Mass_dust = np.zeros(len(gals_dust['Type']))
	Ne_Mass_dust = np.zeros(len(gals_dust['Type']))
	Mg_Mass_dust = np.zeros(len(gals_dust['Type']))
	Si_Mass_dust = np.zeros(len(gals_dust['Type']))
	S_Mass_dust = np.zeros(len(gals_dust['Type']))
	Ca_Mass_dust = np.zeros(len(gals_dust['Type']))
	Fe_Mass_dust = np.zeros(len(gals_dust['Type']))
	O_dustmass_dust = np.zeros(len(gals_dust['Type']))

	
	for i in range(0,len(gals_dust['Type'])):
		H_Mass_dust[i] = gals_dust['ColdGas_elements'][i][0]
		He_Mass_dust[i] = gals_dust['ColdGas_elements'][i][1]
		Cb_Mass_dust[i] = gals_dust['ColdGas_elements'][i][2]
		N_Mass_dust[i] = gals_dust['ColdGas_elements'][i][3]
		O_Mass_dust[i] = gals_dust['ColdGas_elements'][i][4]
		Ne_Mass_dust[i] = gals_dust['ColdGas_elements'][i][5]
		Mg_Mass_dust[i] = gals_dust['ColdGas_elements'][i][6]
		Si_Mass_dust[i] = gals_dust['ColdGas_elements'][i][7]
		S_Mass_dust[i] = gals_dust['ColdGas_elements'][i][8]
		Ca_Mass_dust[i] = gals_dust['ColdGas_elements'][i][9]
		Fe_Mass_dust[i] = gals_dust['ColdGas_elements'][i][10]
		
		O_dustmass_dust[i] = gals_dust['Dust_elements'][i][4]
		
		Stellar_Mass_dust[i] = gals_dust['StellarMass'][i]*1.0E10/0.673
		SFR_dust[i] = gals_dust['Sfr'][i]


	Stellar_Mass_Condition = 8.5
	Dust_Mass_Condition = 0.0 
		
		
		
		
	#---------------------All NEW dust Plots
	
	
	
	
	condition = np.logical_and(np.logical_and(np.logical_and(O_Mass_metals>0,np.log10(SFR_metals)>-2.0),np.log10(SFR_metals)<1.6),np.logical_and(np.log10(Stellar_Mass_metals) > 8.5,np.log10(Stellar_Mass_metals) < 11.5))
	
	
	log_Stellar_Mass_metals = np.log10(Stellar_Mass_metals[condition==1])
	O_metals = O_Mass_metals[condition==1]
	H_metals = H_Mass_metals[condition==1]
	
	condition = np.logical_and(np.logical_and(np.logical_and(O_Mass_dust>0,np.log10(SFR_dust)>-2.0),np.log10(SFR_dust)<1.6),np.logical_and(np.log10(Stellar_Mass_dust) > 8.5,np.log10(Stellar_Mass_dust) < 11.5))
	
	log_Stellar_Mass_dust = np.log10(Stellar_Mass_dust[condition==1])
	O_dust = O_Mass_dust[condition==1]
	H_dust = H_Mass_dust[condition==1]
	Odustdust = O_dustmass_dust[condition==1]
	
	
	log_metallicity_metals = np.log10( (O_metals/H_metals) * (1.0/16.0) )+ 12.0
	log_metallicity_dust   = np.log10( (O_dust/H_dust) * (1.0/16.0) )+ 12.0
	log_metallicity_dust2  = np.log10( ((O_dust+Odustdust)/H_dust) * (1.0/16.0) )+ 12.0
	
	
	SM_bins,Z_bins,Z_std_dev,Z_std_err,count = fit_scatter(log_Stellar_Mass_metals, log_metallicity_metals, ret_n=True, ret_sterr=True, nbins=10)
	SM_bins2,Z_bins2,Z_std_dev2,Z_std_err2,count2 = fit_scatter(log_Stellar_Mass_dust, log_metallicity_dust, ret_n=True, ret_sterr=True, nbins=10)
	SM_bins3,Z_bins3,Z_std_dev3,Z_std_err3,count3 = fit_scatter(log_Stellar_Mass_dust, log_metallicity_dust2, ret_n=True, ret_sterr=True, nbins=10)

	print loop, count

	if(sum(count)>0):	
		plt.xlim([7,12])
		#plt.ylim([7.5,10])
		plt.hexbin(log_Stellar_Mass_metals,log_metallicity_metals,gridsize=500,mincnt=1)
		plt.errorbar(SM_bins,Z_bins,yerr=(Z_std_err),color='r',label='Normal')
		plt.errorbar(SM_bins2,Z_bins2,yerr=(Z_std_err2),color='k',label='Depleted')
		plt.errorbar(SM_bins3,Z_bins3,yerr=(Z_std_err3),color='g',label='Depleted D+M')
		
		plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
		plt.ylabel(r'12 + log$_{10}$(O/H)', fontsize=14,labelpad=0)
		plt.tick_params(axis='both', which='major', labelsize=10)
		plt.tick_params(axis='both', which='minor', labelsize=8)
		#plt.text(6.2,-1,"N = "+str(sum(count)))
		#plt.text(10,2,"z = "+str(loop)+" :New dust")
		plt.legend(loc='center left')
			
		pylab.savefig('./graphs/metals_O_z'+str(loop)+'.png', bbox_inches=0)
		plt.close()
	
		





