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



#### Read in observational data

#### Remy-Ruyer 2014 z=0
#z0_obs_DM, obs_DM_err, obs_SM = np.loadtxt('../observations/z0.txt',unpack=True)
#z0_obs_DM_act_err = obs_DM*(obs_DM_err/100.0)
#Name DustMass DMerror% SMass MHI MHIerror% 12+log(O/H) MH2,mw MH2,Z
Remy_DM, Remy_DMerr, Remy_SM, Remy_MHI, Remy_MHIerr, Remy_Metals, Remy_H2mw, Remy_MH2z = np.loadtxt('../observations/Remy_Ruyer_2014_KINGFISH_z0.txt',unpack=True,comments='#')
Remy_DM_err_actual = Remy_DM*(Remy_DMerr/100.0)


#### Santini 2014 z=0,1,2
# SM DM DM+err DM-err
Santini_SM_z0, Santini_DM_z0, Santini_DMuperr_z0, Santini_DMdownerr_z0 = np.loadtxt('../observations/Santini_2014_z0.txt',unpack=True,comments='#')
Santini_SM_z1, Santini_DM_z1, Santini_DMuperr_z1, Santini_DMdownerr_z1 = np.loadtxt('../observations/Santini_2014_z1.txt',unpack=True,comments='#')
Santini_SM_z2, Santini_DM_z2, Santini_DMuperr_z2, Santini_DMdownerr_z2 = np.loadtxt('../observations/Santini_2014_z2.txt',unpack=True,comments='#')

#### daCunha 2015
#Phot_z z+err z-err SM SM+err SM-err DM DM+err DM-err
daCunha_z,daCunha_zuperr,daCunha_zdownerr,daCunha_SM,daCunha_SMuperr,daCunha_SMdownerr,daCunha_DM,daCunha_DMuperr,daCunha_DMdownerr = np.loadtxt('../observations/daCunha_2015_z_3_4_5.txt',unpack=True,comments='#')


####Mancini2015

Mancini_z, Mancini_SM, Mancini_SMerr, Mancini_DM, Mancini_DMerr = np.loadtxt('../observations/Mancini_2015_z6_z7.txt',unpack=True,comments='#')


#print obs_DM

print "Redshift [Number of galaxies in each mass bin]"	

for loop in range(0,10):
#for loop in range(0,1):

	fin = open('../data/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()

	AGB_Dust_Mass = np.zeros(len(gals['Type']))
	SNII_Dust_Mass = np.zeros(len(gals['Type']))
	SNIa_Dust_Mass = np.zeros(len(gals['Type']))
	Growth_Dust_Mass = np.zeros(len(gals['Type']))
	All_Dust_Mass = np.zeros(len(gals['Type']))


	Stellar_Mass = np.zeros(len(gals['Type']))
	ColdGas = np.zeros(len(gals['Type']))
	SFR = np.zeros(len(gals['Type']))
	Metals = np.zeros(len(gals['Type']))
	Metals_AGB = np.zeros(len(gals['Type']))
	Metals_SNII = np.zeros(len(gals['Type']))
	Metals_SNIA = np.zeros(len(gals['Type']))
	
	H_Mass = np.zeros(len(gals['Type']))
	He_Mass = np.zeros(len(gals['Type']))
	Cb_Mass = np.zeros(len(gals['Type']))
	N_Mass = np.zeros(len(gals['Type']))
	O_Mass = np.zeros(len(gals['Type']))
	Ne_Mass = np.zeros(len(gals['Type']))
	Mg_Mass = np.zeros(len(gals['Type']))
	Si_Mass = np.zeros(len(gals['Type']))
	S_Mass = np.zeros(len(gals['Type']))
	Ca_Mass = np.zeros(len(gals['Type']))
	Fe_Mass = np.zeros(len(gals['Type']))
	
	
	New_Dust_Mass = np.zeros(len(gals['Type']))






	for i in range(0,len(gals['Type'])):
		#print gals['Type'][i]
		AGB_Dust_Mass[i] = (gals['DustMassISM'][i][0]+gals['DustMassISM'][i][1]+gals['DustMassISM'][i][2]+gals['DustMassISM'][i][3])
		SNII_Dust_Mass[i] = (gals['DustMassISM'][i][4]+gals['DustMassISM'][i][5]+gals['DustMassISM'][i][6]+gals['DustMassISM'][i][7])
		SNIa_Dust_Mass[i] = (gals['DustMassISM'][i][8]+gals['DustMassISM'][i][9]+gals['DustMassISM'][i][10]+gals['DustMassISM'][i][11])
		Growth_Dust_Mass[i] = (gals['DustMassISM'][i][12]+gals['DustMassISM'][i][13]+gals['DustMassISM'][i][14]+gals['DustMassISM'][i][15])
		for j in range(0,16):
			All_Dust_Mass[i] += gals['DustMassISM'][i][j]
 		for j in range(0,11):
			New_Dust_Mass[i] += gals['Dust_elements'][i][j]  
			print gals['Dust_elements'][i][j]
		for j in range(2,11):
			Metals[i] += gals['ColdGas_elements'][i][j]

		#print np.log10(Metals[i]), np.log10(New_Dust_Mass[i])
				
				
				
		H_Mass = gals['ColdGas_elements'][0]
		He_Mass = gals['ColdGas_elements'][1]
		Cb_Mass = gals['ColdGas_elements'][2]
		N_Mass = gals['ColdGas_elements'][3]
		O_Mass = gals['ColdGas_elements'][4]
		Ne_Mass = gals['ColdGas_elements'][5]
		Mg_Mass = gals['ColdGas_elements'][6]
		Si_Mass = gals['ColdGas_elements'][7]
		S_Mass = gals['ColdGas_elements'][8]
		Ca_Mass = gals['ColdGas_elements'][9]
		Fe_Mass = gals['ColdGas_elements'][10]
		
		Stellar_Mass[i] = gals['StellarMass'][i]*1.0E10/0.673
		ColdGas[i] = gals['ColdGas'][i]*1.0E10/0.673
		SFR[i] = gals['Sfr'][i]
		
		Metals_AGB  = gals['MetalsColdGas'][i][0]*1.0E10/0.673
		Metals_SNII = gals['MetalsColdGas'][i][1]*1.0E10/0.673
		Metals_SNIA = gals['MetalsColdGas'][i][2]*1.0E10/0.673
		
		#Stellar_Mass_Condition = 1.0E9
		Stellar_Mass_Condition = 0.0
		Dust_Mass_Condition = 0.0 
		
		
	print np.mean(SFR)
		
	#---------------------All NEW dust Plots

	#condition = np.logical_and(New_Dust_Mass>0,Stellar_Mass>Stellar_Mass_Condition)
	condition = np.logical_and(np.logical_and(New_Dust_Mass>0, np.log10(SFR) > -3.0),Stellar_Mass>Stellar_Mass_Condition)
	
	log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
	log_New_Dust_Mass = np.log10(New_Dust_Mass[condition==1])
	
	SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count = fit_scatter(log_Stellar_Mass, log_New_Dust_Mass, ret_n=True, ret_sterr=True, nbins=10)
	print loop, count

	if(sum(count)>0):	
		plt.xlim([6,12])
		plt.ylim([-2,10.2])
		plt.hexbin(log_Stellar_Mass,log_New_Dust_Mass,gridsize=500,mincnt=1, label='Dust All')
		plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='r',label='Dust All')
		plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
		plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
		plt.tick_params(axis='both', which='major', labelsize=10)
		plt.tick_params(axis='both', which='minor', labelsize=8)
		plt.text(6.2,-1,"N = "+str(sum(count)))
		plt.text(10,2,"z = "+str(loop)+" :New dust")
		if(loop == 0):
			plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='g',label='RemyRuyer2014',fmt='o')
			plt.errorbar(Santini_SM_z0, Santini_DM_z0, yerr = (Santini_DMdownerr_z0, Santini_DMuperr_z0), color='r',label='Santini2014',fmt='o')
		if(loop == 1):
			plt.errorbar(Santini_SM_z1, Santini_DM_z1, yerr = (Santini_DMdownerr_z1, Santini_DMuperr_z1), color='g',label='Santini2014',fmt='o')
		if(loop == 2):
			plt.errorbar(Santini_SM_z2, Santini_DM_z2, yerr = (Santini_DMdownerr_z2, Santini_DMuperr_z2), color='g',label='Santini2014',fmt='o')
		if( (loop == 3) or (loop == 4) or (loop == 5) ):
			plt.errorbar(daCunha_SM, daCunha_DM, yerr = (daCunha_DMdownerr, daCunha_DMuperr), color='g',label='daCunha2014',fmt='o')
		if( (loop == 6) or (loop == 7) ):
			plt.errorbar(Mancini_SM, Mancini_DM, yerr = Mancini_DMerr, xerr = Mancini_SMerr, color='g',label='Mancini2015',fmt='o')
		plt.legend(loc='lower right')
			
		pylab.savefig('./graphs/stellar_Newdust_z'+str(loop)+'.png', bbox_inches=0)
		plt.close()
	
	if loop == 0:
		avg_NEW_dust = np.array([np.mean(log_New_Dust_Mass)])
		std_NEW_dust = np.array([np.std(log_New_Dust_Mass)/np.sqrt(len(log_New_Dust_Mass))])

	else:
		avg_NEW_dust = np.append(avg_NEW_dust,np.mean(log_New_Dust_Mass))
		std_NEW_dust = np.append(std_NEW_dust,np.std(log_New_Dust_Mass)/np.sqrt(len(log_New_Dust_Mass)))
		
		

	#---------------------AGB Plots

	condition = np.logical_and(AGB_Dust_Mass>0,Stellar_Mass>Stellar_Mass_Condition)
	
	log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
	log_AGB_Dust_Mass = np.log10(AGB_Dust_Mass[condition==1])

# 	if loop==0:	
 #		for i in range(0,len(log_Stellar_Mass)):
 #			print log_Stellar_Mass[i],log_AGB_Dust_Mass[i]
		
	SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count = fit_scatter(log_Stellar_Mass, log_AGB_Dust_Mass, ret_n=True, ret_sterr=True, nbins=10)
	print loop, count

	if(sum(count)>0):	
		plt.xlim([6,12])
		plt.ylim([-2,10.2])
		plt.hexbin(log_Stellar_Mass,log_AGB_Dust_Mass,gridsize=500,mincnt=1, label='Dust AGB')
		plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='r',label='Dust AGB')
		plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
		plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
		plt.tick_params(axis='both', which='major', labelsize=10)
		plt.tick_params(axis='both', which='minor', labelsize=8)
		#plt.legend(loc='lower right')
		plt.text(6.2,4,"N = "+str(sum(count)))
		plt.text(10,-5,"z = "+str(loop)+" :AGB dust")
	
	
	
		pylab.savefig('./graphs/stellar_AGBdust_z'+str(loop)+'.png', bbox_inches=0)
		plt.close()
	
	if loop == 0:
		redshift = np.array([0])
		avg_AGB_dust = np.array([np.mean(log_AGB_Dust_Mass)])
		std_AGB_dust = np.array([np.std(log_AGB_Dust_Mass)/np.sqrt(len(log_AGB_Dust_Mass))])
		
	else:
		redshift = np.append(redshift,loop)
		avg_AGB_dust = np.append(avg_AGB_dust,np.mean(log_AGB_Dust_Mass))
		std_AGB_dust = np.append(std_AGB_dust,np.std(log_AGB_Dust_Mass)/np.sqrt(len(log_AGB_Dust_Mass)))
	#print redshift
	#print avg_AGB_dust
	#print std_AGB_dust

	#---------------------SNII Plots

	condition = np.logical_and(SNII_Dust_Mass>0,Stellar_Mass>Stellar_Mass_Condition)
	
	log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
	log_SNII_Dust_Mass = np.log10(SNII_Dust_Mass[condition==1])
	
	SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count = fit_scatter(log_Stellar_Mass, log_SNII_Dust_Mass, ret_n=True, ret_sterr=True, nbins=10)
	print loop, count
	
	if(sum(count)>0):	
		plt.xlim([6,12])
		plt.ylim([-2,10.2])
		plt.hexbin(log_Stellar_Mass,log_SNII_Dust_Mass,gridsize=500,mincnt=1, label='Dust SNII')
		plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='r',label='Dust SNII')
		plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
		plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
		plt.tick_params(axis='both', which='major', labelsize=10)
		plt.tick_params(axis='both', which='minor', labelsize=8)
		#plt.legend(loc='lower right')
		plt.text(6.2,4,"N = "+str(sum(count)))
		plt.text(10,-5,"z = "+str(loop)+" :SNII dust")
		pylab.savefig('./graphs/stellar_SNIIdust_z'+str(loop)+'.png', bbox_inches=0)
		plt.close()
	
	if loop == 0:
		avg_SNII_dust = np.array([np.mean(log_SNII_Dust_Mass)])
		std_SNII_dust = np.array([np.std(log_SNII_Dust_Mass)/np.sqrt(len(log_SNII_Dust_Mass))])
	else:
		avg_SNII_dust = np.append(avg_SNII_dust,np.mean(log_SNII_Dust_Mass))
		std_SNII_dust = np.append(std_SNII_dust,np.std(log_SNII_Dust_Mass)/np.sqrt(len(log_SNII_Dust_Mass)))
	
	#---------------------SNIa Plots

	condition = np.logical_and(SNIa_Dust_Mass>0,Stellar_Mass>Stellar_Mass_Condition)
	
	log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
	log_SNIa_Dust_Mass = np.log10(SNIa_Dust_Mass[condition==1])
	
	SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count = fit_scatter(log_Stellar_Mass, log_SNIa_Dust_Mass, ret_n=True, ret_sterr=True, nbins=10)
	print loop, count
	
	if(sum(count)>0):	
		plt.xlim([6,12])
		plt.ylim([-2,10.2])
		plt.hexbin(log_Stellar_Mass,log_SNIa_Dust_Mass,gridsize=500,mincnt=1, label='Dust SNIa')
		plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='r',label='Dust SNIa')
		plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
		plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
		plt.tick_params(axis='both', which='major', labelsize=10)
		plt.tick_params(axis='both', which='minor', labelsize=8)
		#plt.legend(loc='lower right')
		plt.text(6.2,4,"N = "+str(sum(count)))
		plt.text(10,-5,"z = "+str(loop)+" :SNIa dust")
		pylab.savefig('./graphs/stellar_SNIadust_z'+str(loop)+'.png', bbox_inches=0)
		plt.close()
	
	if loop == 0:
		avg_SNIa_dust = np.array([np.mean(log_SNIa_Dust_Mass)])
		std_SNIa_dust = np.array([np.std(log_SNIa_Dust_Mass)/np.sqrt(len(log_SNIa_Dust_Mass))])
	else:
		avg_SNIa_dust = np.append(avg_SNIa_dust,np.mean(log_SNIa_Dust_Mass))
		std_SNIa_dust = np.append(std_SNIa_dust,np.std(log_SNIa_Dust_Mass)/np.sqrt(len(log_SNIa_Dust_Mass)))


	#---------------------Growth Plots

	condition = np.logical_and(Growth_Dust_Mass>0,Stellar_Mass>Stellar_Mass_Condition)
	
	log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
	log_Growth_Dust_Mass = np.log10(Growth_Dust_Mass[condition==1])
	
	SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count = fit_scatter(log_Stellar_Mass, log_Growth_Dust_Mass, ret_n=True, ret_sterr=True, nbins=10)
	print loop, count
	
	if(sum(count)>0):	
		plt.xlim([6,12])
		plt.ylim([-2,10.2])
		plt.hexbin(log_Stellar_Mass,log_Growth_Dust_Mass,gridsize=500,mincnt=1, label='Dust Growth')
		plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='r',label='Dust Growth')
		plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
		plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
		plt.tick_params(axis='both', which='major', labelsize=10)
		plt.tick_params(axis='both', which='minor', labelsize=8)
		#plt.legend(loc='lower right')
		plt.text(6.2,4,"N = "+str(sum(count)))
		plt.text(10,-5,"z = "+str(loop)+" :Growth dust")
		pylab.savefig('./graphs/stellar_GROWdust_z'+str(loop)+'.png', bbox_inches=0)
		plt.close()
	
	if loop == 0:
		avg_GROW_dust = np.array([np.mean(log_Growth_Dust_Mass)])
		std_GROW_dust = np.array([np.std(log_Growth_Dust_Mass)/np.sqrt(len(log_Growth_Dust_Mass))])
		
	else:
		avg_GROW_dust = np.append(avg_GROW_dust,np.mean(log_Growth_Dust_Mass))
		std_GROW_dust = np.append(std_GROW_dust,np.std(log_Growth_Dust_Mass)/np.sqrt(len(log_Growth_Dust_Mass)))
		


	#---------------------All dust Plots

	condition = np.logical_and(All_Dust_Mass>0,Stellar_Mass>Stellar_Mass_Condition)
	
	log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
	log_All_Dust_Mass = np.log10(All_Dust_Mass[condition==1])
	
	SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count = fit_scatter(log_Stellar_Mass, log_All_Dust_Mass, ret_n=True, ret_sterr=True, nbins=10)
	print loop, count

	if(sum(count)>0):	
		plt.xlim([6,12])
		plt.ylim([-2,10.2])
		plt.hexbin(log_Stellar_Mass,log_All_Dust_Mass,gridsize=500,mincnt=1, label='Dust All')
		plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='r',label='Dust All')
		plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
		plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
		plt.tick_params(axis='both', which='major', labelsize=10)
		plt.tick_params(axis='both', which='minor', labelsize=8)
		#plt.legend(loc='lower right')
		plt.text(6.2,4,"N = "+str(sum(count)))
		plt.text(10,-5,"z = "+str(loop)+" :All dust")
		
		pylab.savefig('./graphs/stellar_ALLdust_z'+str(loop)+'.png', bbox_inches=0)
		plt.close()
	
	if loop == 0:
		avg_ALL_dust = np.array([np.mean(log_All_Dust_Mass)])
		std_ALL_dust = np.array([np.std(log_All_Dust_Mass)/np.sqrt(len(log_All_Dust_Mass))])

	else:
		avg_ALL_dust = np.append(avg_ALL_dust,np.mean(log_All_Dust_Mass))
		std_ALL_dust = np.append(std_ALL_dust,np.std(log_All_Dust_Mass)/np.sqrt(len(log_All_Dust_Mass)))




	#---------------------Metal Dust ratios

	condition = np.logical_and(np.logical_and(Metals>0,New_Dust_Mass>0),Stellar_Mass>0)
	
	log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
	log_New_Dust_Mass = np.log10(New_Dust_Mass[condition==1])
	log_Metals = np.log10(Metals[condition==1])
	
	Ratio = log_New_Dust_Mass - log_Metals
	
	SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count = fit_scatter(log_Stellar_Mass, Ratio, ret_n=True, ret_sterr=True, nbins=10)
	print loop, count

	if(sum(count)>0):	
		plt.xlim([6,12])
		plt.ylim([-8,2])
		plt.hexbin(log_Stellar_Mass,Ratio,gridsize=500,mincnt=1, label='Dust All')
		plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='r',label='Dust/Metal ratio')
		plt.xlabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=10)
		plt.ylabel(r'log$_{10}$(Mdust/Mmetals$)', fontsize=14,labelpad=0)
		plt.tick_params(axis='both', which='major', labelsize=10)
		plt.tick_params(axis='both', which='minor', labelsize=8)
		#plt.legend(loc='lower right')
		plt.text(2,-7,"N = "+str(sum(count)))
		plt.text(2,-5,"z = "+str(loop)+" :All dust")
		
		pylab.savefig('./graphs/stellar_dustmetalratio_z'+str(loop)+'.png', bbox_inches=0)
		plt.close()


	#---------------------Gas Dust ratios

	condition = np.logical_and(np.logical_and(ColdGas>0,New_Dust_Mass>0),Stellar_Mass>0)
	
	log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
	log_New_Dust_Mass = np.log10(New_Dust_Mass[condition==1])
	log_ColdGas = np.log10(ColdGas[condition==1])
	
	Ratio = log_New_Dust_Mass - log_ColdGas
	
	SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count = fit_scatter(log_Stellar_Mass, Ratio, ret_n=True, ret_sterr=True, nbins=10)
	print loop, count

	if(sum(count)>0):	
		plt.xlim([6,12])
		plt.ylim([-8,2])
		plt.hexbin(log_Stellar_Mass,Ratio,gridsize=500,mincnt=1, label='Dust All')
		plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='r',label='Dust/Metal ratio')
		plt.xlabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=10)
		plt.ylabel(r'log$_{10}$(Mdust/Mcoldgas$)', fontsize=14,labelpad=0)
		plt.tick_params(axis='both', which='major', labelsize=10)
		plt.tick_params(axis='both', which='minor', labelsize=8)
		#plt.legend(loc='lower right')
		plt.text(2,-7,"N = "+str(sum(count)))
		plt.text(2,-5,"z = "+str(loop)+" :All dust")
		
		pylab.savefig('./graphs/stellar_dustgasratio_z'+str(loop)+'.png', bbox_inches=0)
		plt.close()


	#---------------------Dust mass functions

	
	log_New_Dust_Mass = np.log10(New_Dust_Mass[New_Dust_Mass>0.0])
	
	hist, bin_edges = np.histogram(log_New_Dust_Mass)	
	bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
	
	print hist
	print bin_edges
	
	plt.xlim([6,12])
	#plt.ylim([-8,2])
	
 	volume = (480.279/0.673)**3.0
# 	volume = (0.938044921/0.673)**3.0
	binsize = 1.35135751
	
	plt.plot(bin_centers,hist/(volume*binsize),color='k',label='Dust All')
	
	plt.xlabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=10)
	#plt.ylabel(r'log$_{10}$(Mdust/Mcoldgas$)', fontsize=14,labelpad=0)
	plt.tick_params(axis='both', which='major', labelsize=10)
	plt.tick_params(axis='both', which='minor', labelsize=8)
	#plt.legend(loc='lower right')
	#plt.text(2,-7,"N = "+str(sum(count)))
	#plt.text(2,-5,"z = "+str(loop)+" :All dust")
	
	pylab.savefig('./graphs/dustmass_function_z'+str(loop)+'.png', bbox_inches=0)
	plt.close()



#---------------------Redshift vs Dust mass plot
# print redshift
# print avg_AGB_dust
# print std_AGB_dust
# print avg_SNII_dust
# print std_SNII_dust
# print avg_GROW_dust
# print std_GROW_dust
# print avg_ALL_dust
# print std_ALL_dust


plt.xlabel(r'redshift', fontsize=14,labelpad=10)
plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
plt.xlim(-1,12)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.errorbar(redshift,avg_AGB_dust ,std_AGB_dust, color='r',label='AGB')
plt.errorbar(redshift,avg_SNII_dust,std_SNII_dust, color='b',label='SNII')
plt.errorbar(redshift,avg_GROW_dust,std_GROW_dust, color='g',label='Grown')
plt.errorbar(redshift,avg_SNIa_dust,std_SNIa_dust, color='y',label='SNIA')
plt.errorbar(redshift,avg_ALL_dust,std_ALL_dust, color='k',label='ALL')
#plt.errorbar(redshift,avg_NEW_dust,std_NEW_dust, color='yellow',label='New')

plt.legend(loc='upper right')
pylab.savefig('./graphs/dust_redshift.png', bbox_inches=0)




# 0		dust.AGB.SiC    
# 1		dust.AGB.Sil    
# 2		dust.AGB.Cb     
# 3		dust.AGB.Fe    
 
# 4		dust.SNII.SiC   
# 5		dust.SNII.Sil   
# 6		dust.SNII.Cb    
# 7		dust.SNII.Fe    

# 8		dust.SNIa.SiC   
# 9		dust.SNIa.Sil   
# 10	dust.SNIa.Cb    
# 11	dust.SNIa.Fe    

# 12	dust.Growth.SiC 
# 13	dust.Growth.Sil 
# 14	dust.Growth.Cb  
# 15	dust.Growth.Fe  

















