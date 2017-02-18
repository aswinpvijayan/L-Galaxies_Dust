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





#### Santini 2014 z=0,1,2
# SM DM DM+err DM-err
# Santini_SM_z0, Santini_DM_z0, Santini_DMuperr_z0, Santini_DMdownerr_z0 = np.loadtxt('../observations/Santini_2014_z0.txt',unpack=True,comments='#')


#for loop in range(0,10):
for loop in range(4,5):

	fin = open('../data/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()

	#Type = np.zeros(len(gals['Type']))
	Type = gals['Type']
	#print Type
	
	UV_mag = gals['Mag'][27]
	UV_magdust = gals['MagDust'][27]
	
	print UV_mag
	print UV_magdust
	
	print UV_mag - UV_magdust
	

# def luminosity_function(G_MR, ThisRedshiftList, Volume):
           
	FUV_lim=[-23.0,-16.0]
	phi_lim=[-7.,2.]

	BoxSize = 96.0558/512.0 #1/512 MRII
	BoxSize = 480.279/512.0 #1/512 MR
	
	Volume = BoxSize**3.0
	bin=0.5
	Nbins=int((FUV_lim[1]-FUV_lim[0])/bin)
	bin_arr=np.arange(FUV_lim[0],FUV_lim[1]+bin,bin)


	hist=np.histogram(UV_mag, bins=bin_arr, range=(FUV_lim[0],FUV_lim[1]))
	hist2=np.histogram(UV_magdust, bins=bin_arr, range=(FUV_lim[0],FUV_lim[1]))
	
	print hist[1]
	print np.log10(hist[0]/(Volume*bin))
	
	fig=plt.figure()
	ax=fig.add_subplot(111)
	ax.plot(hist[1][0:len(hist[1][:])-1]+bin/2.,np.log10(hist[0][:]/(Volume*bin)),color='red', linewidth=2)
	ax.plot(hist2[1][0:len(hist2[1][:])-1]+bin/2.,np.log10(hist2[0][:]/(Volume*bin)),color='blue', linewidth=2)
	ax.set_xlabel('$\mathrm{Muv}$')
	ax.set_ylabel('$\mathrm{log_{10}}(\phi [h^3 \mathrm{Mpc^{-3}} \mathrm{mag^{-1}}])$')
	ax.set_xlim(FUV_lim[0],FUV_lim[1])
	ax.set_ylim(phi_lim[0],phi_lim[1])
	#ax.xaxis.set_major_locator(MultipleLocator(1))
	#ax.yaxis.set_minor_locator(MultipleLocator(0.5))
	#ax.text(-22.,-1.,'z='+str(z_in),fontsize=14)
	fig.subplots_adjust(bottom=0.15)

	fname = 'LF_z'+str(loop)+'.png'
	fig.savefig(fname)
	fig.clf()






# 	for i in range(0,len(gals['Type'])):
				
				
				
# 		H_Mass = gals['ColdGas_elements'][0]
# 		He_Mass = gals['ColdGas_elements'][1]
# 		Cb_Mass = gals['ColdGas_elements'][2]
# 		N_Mass = gals['ColdGas_elements'][3]
# 		O_Mass = gals['ColdGas_elements'][4]
# 		Ne_Mass = gals['ColdGas_elements'][5]
# 		Mg_Mass = gals['ColdGas_elements'][6]
# 		Si_Mass = gals['ColdGas_elements'][7]
# 		S_Mass = gals['ColdGas_elements'][8]
# 		Ca_Mass = gals['ColdGas_elements'][9]
# 		Fe_Mass = gals['ColdGas_elements'][10]
# 		
# 		Stellar_Mass[i] = gals['StellarMass'][i]*1.0E10/0.673
# 		ColdGas[i] = gals['ColdGas'][i]*1.0E10/0.673
# 		SFR[i] = gals['Sfr'][i]
# 		Type[i] = gals['Type'][i]
# 		sSFR[i] = SFR[i] / Stellar_Mass[i]
# 		#print sSFR[i]
# 		
# 		Metals_AGB  = gals['MetalsColdGas'][i][0]*1.0E10/0.673
# 		Metals_SNII = gals['MetalsColdGas'][i][1]*1.0E10/0.673
# 		Metals_SNIA = gals['MetalsColdGas'][i][2]*1.0E10/0.673
# 		

# 		
# 	#---------------------All NEW dust Plots
# 
# 	#condition = np.logical_and(New_Dust_Mass>0,Stellar_Mass>Stellar_Mass_Condition)
# 	#condition = np.logical_and(np.logical_and(New_Dust_Mass>0, np.log10(SFR) > -3.0),Stellar_Mass>Stellar_Mass_Condition)
# 	condition = np.logical_and(sSFR > 0.0345E-9, np.logical_and(New_Dust_Mass > 0, Type == 0))
# 	
# 	log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
# 	log_New_Dust_Mass = np.log10(New_Dust_Mass[condition==1])
# 	
# 	SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count = fit_scatter(log_Stellar_Mass, log_New_Dust_Mass, ret_n=True, ret_sterr=True, nbins=10)
# 	print loop, count
# 
# 	if(sum(count)>0):	
# 		plt.xlim([6,12])
# 		plt.ylim([-2,10.2])
# 		plt.hexbin(log_Stellar_Mass,log_New_Dust_Mass,gridsize=500,mincnt=1, label='Dust All')
# 		plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='r',label='Dust All')
# 		plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
# 		plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
# 		plt.tick_params(axis='both', which='major', labelsize=10)
# 		plt.tick_params(axis='both', which='minor', labelsize=8)
# 		plt.text(6.2,-1,"N = "+str(sum(count)))
# 		plt.text(10,2,"z = "+str(loop)+" :New dust")
# 		if(loop == 0):
# 			plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='g',label='RemyRuyer2014',fmt='o')
# 			plt.errorbar(Santini_SM_z0, Santini_DM_z0, yerr = (Santini_DMdownerr_z0, Santini_DMuperr_z0), color='r',label='Santini2014',fmt='o')
# 		if(loop == 1):
# 			plt.errorbar(Santini_SM_z1, Santini_DM_z1, yerr = (Santini_DMdownerr_z1, Santini_DMuperr_z1), color='g',label='Santini2014',fmt='o')
# 		if(loop == 2):
# 			plt.errorbar(Santini_SM_z2, Santini_DM_z2, yerr = (Santini_DMdownerr_z2, Santini_DMuperr_z2), color='g',label='Santini2014',fmt='o')
# 		if( (loop == 3) or (loop == 4) or (loop == 5) ):
# 			plt.errorbar(daCunha_SM, daCunha_DM, yerr = (daCunha_DMdownerr, daCunha_DMuperr), color='g',label='daCunha2014',fmt='o')
# 		if( (loop == 6) or (loop == 7) ):
# 			plt.errorbar(Mancini_SM, Mancini_DM, yerr = Mancini_DMerr, xerr = Mancini_SMerr, color='g',label='Mancini2015',fmt='o')
# 		plt.legend(loc='lower right')
# 			
# 		pylab.savefig('./graphs/stellar_Newdust_z'+str(loop)+'.png', bbox_inches=0)
# 		plt.close()
# 	
# 	if loop == 0:
# 		avg_NEW_dust = np.array([np.mean(log_New_Dust_Mass)])
# 		std_NEW_dust = np.array([np.std(log_New_Dust_Mass)/np.sqrt(len(log_New_Dust_Mass))])
# 
# 	else:
# 		avg_NEW_dust = np.append(avg_NEW_dust,np.mean(log_New_Dust_Mass))
# 		std_NEW_dust = np.append(std_NEW_dust,np.std(log_New_Dust_Mass)/np.sqrt(len(log_New_Dust_Mass)))
# 		
# 		
