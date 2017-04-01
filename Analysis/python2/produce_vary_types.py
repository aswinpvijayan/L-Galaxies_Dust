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


# 
# #####Colour scheme
# tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
#              (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
#              (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
#              (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
#              (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
#   
# # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
# for i in range(len(tableau20)):    
#     r, g, b = tableau20[i]    
#     tableau20[i] = (r / 255., g / 255., b / 255.)    
 

#### Read in observational data

#### Remy-Ruyer 2014 z=0
#z0_obs_DM, obs_DM_err, obs_SM = np.loadtxt('../observations/z0.txt',unpack=True)
#z0_obs_DM_act_err = obs_DM*(obs_DM_err/100.0)
#Name DustMass DMerror% SMass MHI MHIerror% 12+log(O/H) MH2,mw MH2,Z
Remy_DM, Remy_DMerr, Remy_SM, Remy_MHI, Remy_MHIerr, Remy_Metals, Remy_H2mw, Remy_MH2z = np.loadtxt('../observations/Remy_Ruyer_2014_KINGFISH_z0.txt',unpack=True,comments='#')
Remy_DM_err_actual = Remy_DM*(Remy_DMerr/100.0)
Remy_Dust_Gas_Ratio = Remy_DM / (Remy_MHI + Remy_H2mw)
Remy_Dust_Metal_Ratio = Remy_DM / (((10**(Remy_Metals - 8.69))*0.0134)*(Remy_MHI + Remy_H2mw))

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


###Ciesla2014
Ciesla_ID1, Ciesla_DM, Ciesla_DMerr, Ciesla_ID2, Ciesla_SM = np.loadtxt('../observations/Ciesla_2014_z0.txt',unpack=True,comments='#')


# for i in range(len(Ciesla_ID1)):
# 	if Ciesla_ID1[i] != Ciesla_ID2[i]:
# 		print Ciesla_ID1[i], Ciesla_ID2[i]

#### Wiseman2016
Wiseman_z_z2, Wiseman_SM_z2, Wiseman_SMuperr_z2, Wiseman_SMdownerr_z2, Wiseman_SFR_z2, Wiseman_SFRerr_z2, Wiseman_Metals_z2, Wiseman_Metalserr_z2, Wiseman_DTM_z2, Wiseman_DTMerr_z2 = np.loadtxt('../observations/wiseman2016_z2.txt',unpack=True,comments='#')
Wiseman_z_z3, Wiseman_SM_z3, Wiseman_SMuperr_z3, Wiseman_SMdownerr_z3, Wiseman_SFR_z3, Wiseman_SFRerr_z3, Wiseman_Metals_z3, Wiseman_Metalserr_z3, Wiseman_DTM_z3, Wiseman_DTMerr_z3 = np.loadtxt('../observations/wiseman2016_z3.txt',unpack=True,comments='#')
Wiseman_z_z4, Wiseman_SM_z4, Wiseman_SMuperr_z4, Wiseman_SMdownerr_z4, Wiseman_SFR_z4, Wiseman_SFRerr_z4, Wiseman_Metals_z4, Wiseman_Metalserr_z4, Wiseman_DTM_z4, Wiseman_DTMerr_z4 = np.loadtxt('../observations/wiseman2016_z4.txt',unpack=True,comments='#')


SM_bins=np.zeros((6,10,10))
Dust_bins=np.zeros((6,10,10))
Dust_std_dev=np.zeros((6,10,10))
Dust_std_err=np.zeros((6,10,10))
count=np.zeros((6,10,10))




print "Redshift [Number of galaxies in each mass bin]"	
for loop in range(0,10):
#for loop in range(0,1):

	#Read in L-Galaxies data----------------------------
	#THE MODEL
	#tacc = 0.2Myr
	#tdes = 500Myr
	fin = open('../data/MR_tacc02_tdes500/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()

	Type_1 = np.zeros(len(gals['Type']))
	Stellar_Mass_1 = np.zeros(len(gals['Type']))
	ColdGas_1 = np.zeros(len(gals['Type']))
	
	Dust_Mass_1 = np.zeros(len(gals['Type']))
	for i in range(0,len(gals['Type'])):
 		for j in range(0,11):
			Dust_Mass_1[i] += gals['Dust_elements'][i][j]  		
	Stellar_Mass_1 = gals['StellarMass']*1.0E10/0.673
	ColdGas_1 = gals['ColdGas']*1.0E10/0.673
	SFR_1 = gals['Sfr']
	Type_1 = gals['Type']
	sSFR_1 = SFR_1 / Stellar_Mass_1
		
	
	#Read in L-Galaxies data----------------------------
	#tacc = 2Myr
	#tdes = 500Myr
	fin = open('../data/MR_AGB_only/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()

	Type_2 = np.zeros(len(gals['Type']))
	Stellar_Mass_2 = np.zeros(len(gals['Type']))
	ColdGas_2 = np.zeros(len(gals['Type']))
	
	Dust_Mass_2 = np.zeros(len(gals['Type']))
	for i in range(0,len(gals['Type'])):
 		for j in range(0,11):
			Dust_Mass_2[i] += gals['Dust_elements'][i][j]  		
	Stellar_Mass_2 = gals['StellarMass']*1.0E10/0.673
	ColdGas_2 = gals['ColdGas']*1.0E10/0.673
	SFR_2 = gals['Sfr']
	Type_2 = gals['Type']
	sSFR_2 = SFR_2 / Stellar_Mass_2
	
	
	#Read in L-Galaxies data----------------------------
	#tacc = 20Myr
	#tdes = 500Myr
	fin = open('../data/MR_SNII_only/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()

	Type_3 = np.zeros(len(gals['Type']))
	Stellar_Mass_3 = np.zeros(len(gals['Type']))
	ColdGas_3 = np.zeros(len(gals['Type']))
	
	Dust_Mass_3 = np.zeros(len(gals['Type']))
	for i in range(0,len(gals['Type'])):
 		for j in range(0,11):
			Dust_Mass_3[i] += gals['Dust_elements'][i][j]  		
	Stellar_Mass_3 = gals['StellarMass']*1.0E10/0.673
	ColdGas_3 = gals['ColdGas']*1.0E10/0.673
	SFR_3 = gals['Sfr']
	Type_3 = gals['Type']
	sSFR_3 = SFR_3 / Stellar_Mass_3
	
	#Read in L-Galaxies data----------------------------
	#tacc = 0.2Myr
	#tdes = 50Myr
	fin = open('../data/MR_SNIA_only/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()

	Type_4 = np.zeros(len(gals['Type']))
	Stellar_Mass_4 = np.zeros(len(gals['Type']))
	ColdGas_4 = np.zeros(len(gals['Type']))
	
	Dust_Mass_4 = np.zeros(len(gals['Type']))
	for i in range(0,len(gals['Type'])):
 		for j in range(0,11):
			Dust_Mass_4[i] += gals['Dust_elements'][i][j]  		
	Stellar_Mass_4 = gals['StellarMass']*1.0E10/0.673
	ColdGas_4 = gals['ColdGas']*1.0E10/0.673
	SFR_4 = gals['Sfr']
	Type_4 = gals['Type']
	sSFR_4 = SFR_4 / Stellar_Mass_4
	
	
	#Read in L-Galaxies data----------------------------
	#tacc = 0.2Myr
	#tdes = 5000Myr
	fin = open('../data/MR_GROW_only/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()

	Type_5 = np.zeros(len(gals['Type']))
	Stellar_Mass_5 = np.zeros(len(gals['Type']))
	ColdGas_5 = np.zeros(len(gals['Type']))
	
	Dust_Mass_5 = np.zeros(len(gals['Type']))
	for i in range(0,len(gals['Type'])):
 		for j in range(0,11):
			Dust_Mass_5[i] += gals['Dust_elements'][i][j]  		
	Stellar_Mass_5 = gals['StellarMass']*1.0E10/0.673
	ColdGas_5 = gals['ColdGas']*1.0E10/0.673
	SFR_5 = gals['Sfr']
	Type_5 = gals['Type']
	sSFR_5 = SFR_5 / Stellar_Mass_5
	
	#Read in L-Galaxies data----------------------------
	#tacc = 0.2Myr
	#tdes = 5000Myr
	fin = open('../data/MR_ALL_NO_DEST/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()

	Type_6 = np.zeros(len(gals['Type']))
	Stellar_Mass_6 = np.zeros(len(gals['Type']))
	ColdGas_6 = np.zeros(len(gals['Type']))
	
	Dust_Mass_6 = np.zeros(len(gals['Type']))
	for i in range(0,len(gals['Type'])):
 		for j in range(0,11):
			Dust_Mass_6[i] += gals['Dust_elements'][i][j]  		
	Stellar_Mass_6 = gals['StellarMass']*1.0E10/0.673
	ColdGas_6 = gals['ColdGas']*1.0E10/0.673
	SFR_6 = gals['Sfr']
	Type_6 = gals['Type']
	sSFR_6 = SFR_6 / Stellar_Mass_6





	#---------------------SM vs. DM for varied tacc_0 rates

	fig = plt.figure(figsize=(9,13.5))

	#AGB
	condition = np.logical_and(sSFR_2 > 0.0345E-9, np.logical_and(Dust_Mass_2 > 0, Type_2 == 0))
	log_Stellar_Mass_2 = np.log10(Stellar_Mass_2[condition==1])
	log_New_Dust_Mass_2 = np.log10(Dust_Mass_2[condition==1])
	SM_bins[1][loop],Dust_bins[1][loop],Dust_std_dev[1][loop],Dust_std_err[1][loop],count[1][loop] = fit_scatter(log_Stellar_Mass_2, log_New_Dust_Mass_2, ret_n=True, ret_sterr=True, nbins=10)
	plt.subplot(3,2,1)
	plt.xlim([8,12])
	plt.ylim([0,10.2])
	plt.hexbin(log_Stellar_Mass_2,log_New_Dust_Mass_2,gridsize=500,mincnt=1, label='L-Galaxies 2Dhist')
	plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
	plt.ylabel(r'log$_{10}$(M$_{\rm{d}}$/M$_{\odot}$)', fontsize=18)
	plt.tick_params(axis='both', which='major', labelsize=12,width=2)
	plt.tick_params(axis='both', which='minor', labelsize=12,width=2)
	#plt.text(7,1,"N = "+str(sum(count))+"\nz = "+str(loop))
	plt.text(9.6,1,"z = "+str(loop)+"\nAGB",fontsize=16)
	if(loop == 0):
		plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='g',label='RemyRuyer2014',fmt='o')
		plt.errorbar(np.log10(Bourne_MEDMSTAR), np.log10(Bourne_MEDMDUST), yerr = np.log10(Bourne_MEDMDUSTERR/Bourne_MEDMDUST) , color='orange',label='Bourne2012',fmt='o')
		plt.errorbar(Ciesla_SM, Ciesla_DM, yerr = Ciesla_DMerr , color='r',label='Ciesla2014',fmt='o')
		plt.errorbar(Santini_SM_z0, Santini_DM_z0, yerr = (Santini_DMdownerr_z0, Santini_DMuperr_z0), color='b',label='Santini2014',fmt='o')
	if(loop == 1):
		plt.errorbar(Santini_SM_z1, Santini_DM_z1, yerr = (Santini_DMdownerr_z1, Santini_DMuperr_z1), color='g',label='Santini2014',fmt='o')
	if(loop == 2):
		plt.errorbar(Santini_SM_z2, Santini_DM_z2, yerr = (Santini_DMdownerr_z2, Santini_DMuperr_z2), color='g',label='Santini2014',fmt='o')
		plt.errorbar(daCunha_SM_z2, daCunha_DM_z2, yerr = (daCunha_DMdownerr_z2, daCunha_DMuperr_z2), color='b',label='daCunha2014',fmt='o')
	if(loop == 3):
		plt.errorbar(daCunha_SM_z3, daCunha_DM_z3, yerr = (daCunha_DMdownerr_z3, daCunha_DMuperr_z3), color='g',label='daCunha2014',fmt='o')
	if(loop == 4):
		plt.errorbar(daCunha_SM_z4, daCunha_DM_z4, yerr = (daCunha_DMdownerr_z4, daCunha_DMuperr_z4), color='g',label='daCunha2014',fmt='o')
	if(loop == 5):
		plt.errorbar(daCunha_SM_z5, daCunha_DM_z5, yerr = (daCunha_DMdownerr_z5, daCunha_DMuperr_z5), color='g',label='daCunha2014',fmt='o')
	if(loop == 6):
		plt.errorbar(daCunha_SM_z6, daCunha_DM_z6, yerr = (daCunha_DMdownerr_z6, daCunha_DMuperr_z6), color='g',label='daCunha2014',fmt='o')
	if( (loop == 6) or (loop == 7) ):
		plt.errorbar(Mancini_SM, Mancini_DM, yerr = Mancini_DMerr, xerr = Mancini_SMerr, color='g',label='Mancini2015',fmt='o')
	plt.errorbar(SM_bins[1][loop],Dust_bins[1][loop],yerr=(Dust_std_err[0][loop]),color='k',label='tacc=0.2Myr tdes=tdes',linewidth=2)

	#SNII
	condition = np.logical_and(sSFR_3 > 0.0345E-9, np.logical_and(Dust_Mass_3 > 0, Type_3 == 0))
	log_Stellar_Mass_3 = np.log10(Stellar_Mass_3[condition==1])
	log_New_Dust_Mass_3 = np.log10(Dust_Mass_3[condition==1])
	SM_bins[2][loop],Dust_bins[2][loop],Dust_std_dev[2][loop],Dust_std_err[2][loop],count[2][loop] = fit_scatter(log_Stellar_Mass_3, log_New_Dust_Mass_3, ret_n=True, ret_sterr=True, nbins=10)
	plt.subplot(3,2,2)
	plt.xlim([8,12])
	plt.ylim([0,10.2])
	plt.hexbin(log_Stellar_Mass_3,log_New_Dust_Mass_3,gridsize=500,mincnt=1, label='L-Galaxies 2Dhist')
	plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
	plt.ylabel(r'log$_{10}$(M$_{\rm{d}}$/M$_{\odot}$)', fontsize=18)
	plt.tick_params(axis='both', which='major', labelsize=12,width=2)
	plt.tick_params(axis='both', which='minor', labelsize=12,width=2)
	#plt.text(7,1,"N = "+str(sum(count))+"\nz = "+str(loop))
	plt.text(9.6,1,"z = "+str(loop)+"\nSNII",fontsize=16)
	if(loop == 0):
		plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='g',label='RemyRuyer2014',fmt='o')
		plt.errorbar(np.log10(Bourne_MEDMSTAR), np.log10(Bourne_MEDMDUST), yerr = np.log10(Bourne_MEDMDUSTERR/Bourne_MEDMDUST) , color='orange',label='Bourne2012',fmt='o')
		plt.errorbar(Ciesla_SM, Ciesla_DM, yerr = Ciesla_DMerr , color='r',label='Ciesla2014',fmt='o')
		plt.errorbar(Santini_SM_z0, Santini_DM_z0, yerr = (Santini_DMdownerr_z0, Santini_DMuperr_z0), color='b',label='Santini2014',fmt='o')
	if(loop == 1):
		plt.errorbar(Santini_SM_z1, Santini_DM_z1, yerr = (Santini_DMdownerr_z1, Santini_DMuperr_z1), color='g',label='Santini2014',fmt='o')
	if(loop == 2):
		plt.errorbar(Santini_SM_z2, Santini_DM_z2, yerr = (Santini_DMdownerr_z2, Santini_DMuperr_z2), color='g',label='Santini2014',fmt='o')
		plt.errorbar(daCunha_SM_z2, daCunha_DM_z2, yerr = (daCunha_DMdownerr_z2, daCunha_DMuperr_z2), color='b',label='daCunha2014',fmt='o')
	if(loop == 3):
		plt.errorbar(daCunha_SM_z3, daCunha_DM_z3, yerr = (daCunha_DMdownerr_z3, daCunha_DMuperr_z3), color='g',label='daCunha2014',fmt='o')
	if(loop == 4):
		plt.errorbar(daCunha_SM_z4, daCunha_DM_z4, yerr = (daCunha_DMdownerr_z4, daCunha_DMuperr_z4), color='g',label='daCunha2014',fmt='o')
	if(loop == 5):
		plt.errorbar(daCunha_SM_z5, daCunha_DM_z5, yerr = (daCunha_DMdownerr_z5, daCunha_DMuperr_z5), color='g',label='daCunha2014',fmt='o')
	if(loop == 6):
		plt.errorbar(daCunha_SM_z6, daCunha_DM_z6, yerr = (daCunha_DMdownerr_z6, daCunha_DMuperr_z6), color='g',label='daCunha2014',fmt='o')
	if( (loop == 6) or (loop == 7) ):
		plt.errorbar(Mancini_SM, Mancini_DM, yerr = Mancini_DMerr, xerr = Mancini_SMerr, color='g',label='Mancini2015',fmt='o')
	plt.errorbar(SM_bins[2][loop],Dust_bins[2][loop],yerr=(Dust_std_err[0][loop]),color='k',label='tacc=0.2Myr tdes=tdes',linewidth=2)

	#SNIA
	condition = np.logical_and(sSFR_4 > 0.0345E-9, np.logical_and(Dust_Mass_4 > 0, Type_4 == 0))
	log_Stellar_Mass_4 = np.log10(Stellar_Mass_4[condition==1])
	log_New_Dust_Mass_4 = np.log10(Dust_Mass_4[condition==1])
	SM_bins[3][loop],Dust_bins[3][loop],Dust_std_dev[3][loop],Dust_std_err[3][loop],count[3][loop] = fit_scatter(log_Stellar_Mass_4, log_New_Dust_Mass_4, ret_n=True, ret_sterr=True, nbins=10)
	plt.subplot(3,2,3)
	plt.xlim([8,12])
	plt.ylim([0,10.2])
	plt.hexbin(log_Stellar_Mass_4,log_New_Dust_Mass_4,gridsize=500,mincnt=1, label='L-Galaxies 2Dhist')
	plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
	plt.ylabel(r'log$_{10}$(M$_{\rm{d}}$/M$_{\odot}$)', fontsize=18)
	plt.tick_params(axis='both', which='major', labelsize=12,width=2)
	plt.tick_params(axis='both', which='minor', labelsize=12,width=2)
	#plt.text(7,1,"N = "+str(sum(count))+"\nz = "+str(loop))
	plt.text(9.6,1,"z = "+str(loop)+"\nSNIA",fontsize=16)
	if(loop == 0):
		plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='g',label='RemyRuyer2014',fmt='o')
		plt.errorbar(np.log10(Bourne_MEDMSTAR), np.log10(Bourne_MEDMDUST), yerr = np.log10(Bourne_MEDMDUSTERR/Bourne_MEDMDUST) , color='orange',label='Bourne2012',fmt='o')
		plt.errorbar(Ciesla_SM, Ciesla_DM, yerr = Ciesla_DMerr , color='r',label='Ciesla2014',fmt='o')
		plt.errorbar(Santini_SM_z0, Santini_DM_z0, yerr = (Santini_DMdownerr_z0, Santini_DMuperr_z0), color='b',label='Santini2014',fmt='o')
	if(loop == 1):
		plt.errorbar(Santini_SM_z1, Santini_DM_z1, yerr = (Santini_DMdownerr_z1, Santini_DMuperr_z1), color='g',label='Santini2014',fmt='o')
	if(loop == 2):
		plt.errorbar(Santini_SM_z2, Santini_DM_z2, yerr = (Santini_DMdownerr_z2, Santini_DMuperr_z2), color='g',label='Santini2014',fmt='o')
		plt.errorbar(daCunha_SM_z2, daCunha_DM_z2, yerr = (daCunha_DMdownerr_z2, daCunha_DMuperr_z2), color='b',label='daCunha2014',fmt='o')
	if(loop == 3):
		plt.errorbar(daCunha_SM_z3, daCunha_DM_z3, yerr = (daCunha_DMdownerr_z3, daCunha_DMuperr_z3), color='g',label='daCunha2014',fmt='o')
	if(loop == 4):
		plt.errorbar(daCunha_SM_z4, daCunha_DM_z4, yerr = (daCunha_DMdownerr_z4, daCunha_DMuperr_z4), color='g',label='daCunha2014',fmt='o')
	if(loop == 5):
		plt.errorbar(daCunha_SM_z5, daCunha_DM_z5, yerr = (daCunha_DMdownerr_z5, daCunha_DMuperr_z5), color='g',label='daCunha2014',fmt='o')
	if(loop == 6):
		plt.errorbar(daCunha_SM_z6, daCunha_DM_z6, yerr = (daCunha_DMdownerr_z6, daCunha_DMuperr_z6), color='g',label='daCunha2014',fmt='o')
	if( (loop == 6) or (loop == 7) ):
		plt.errorbar(Mancini_SM, Mancini_DM, yerr = Mancini_DMerr, xerr = Mancini_SMerr, color='g',label='Mancini2015',fmt='o')
	plt.errorbar(SM_bins[3][loop],Dust_bins[3][loop],yerr=(Dust_std_err[0][loop]),color='k',label='tacc=0.2Myr tdes=tdes',linewidth=2)

	#GROW
	condition = np.logical_and(sSFR_5 > 0.0345E-9, np.logical_and(Dust_Mass_5 > 10000, Type_5 == 0))
	log_Stellar_Mass_5 = np.log10(Stellar_Mass_5[condition==1])
	log_New_Dust_Mass_5 = np.log10(Dust_Mass_5[condition==1])
	SM_bins[4][loop],Dust_bins[4][loop],Dust_std_dev[4][loop],Dust_std_err[4][loop],count[4][loop] = fit_scatter(log_Stellar_Mass_5, log_New_Dust_Mass_5, ret_n=True, ret_sterr=True, nbins=10)
	plt.subplot(3,2,4)
	plt.xlim([8,12])
	plt.ylim([0,10.2])
	plt.hexbin(log_Stellar_Mass_5,log_New_Dust_Mass_5,gridsize=500,mincnt=1, label='L-Galaxies 2Dhist')
	plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
	plt.ylabel(r'log$_{10}$(M$_{\rm{d}}$/M$_{\odot}$)', fontsize=18)
	plt.tick_params(axis='both', which='major', labelsize=12,width=2)
	plt.tick_params(axis='both', which='minor', labelsize=12,width=2)
	#plt.text(7,1,"N = "+str(sum(count))+"\nz = "+str(loop))
	plt.text(9.6,1,"z = "+str(loop)+"\nMC Growth",fontsize=16)
	if(loop == 0):
		plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='g',label='RemyRuyer2014',fmt='o')
		plt.errorbar(np.log10(Bourne_MEDMSTAR), np.log10(Bourne_MEDMDUST), yerr = np.log10(Bourne_MEDMDUSTERR/Bourne_MEDMDUST) , color='orange',label='Bourne2012',fmt='o')
		plt.errorbar(Ciesla_SM, Ciesla_DM, yerr = Ciesla_DMerr , color='r',label='Ciesla2014',fmt='o')
		plt.errorbar(Santini_SM_z0, Santini_DM_z0, yerr = (Santini_DMdownerr_z0, Santini_DMuperr_z0), color='b',label='Santini2014',fmt='o')
	if(loop == 1):
		plt.errorbar(Santini_SM_z1, Santini_DM_z1, yerr = (Santini_DMdownerr_z1, Santini_DMuperr_z1), color='g',label='Santini2014',fmt='o')
	if(loop == 2):
		plt.errorbar(Santini_SM_z2, Santini_DM_z2, yerr = (Santini_DMdownerr_z2, Santini_DMuperr_z2), color='g',label='Santini2014',fmt='o')
		plt.errorbar(daCunha_SM_z2, daCunha_DM_z2, yerr = (daCunha_DMdownerr_z2, daCunha_DMuperr_z2), color='b',label='daCunha2014',fmt='o')
	if(loop == 3):
		plt.errorbar(daCunha_SM_z3, daCunha_DM_z3, yerr = (daCunha_DMdownerr_z3, daCunha_DMuperr_z3), color='g',label='daCunha2014',fmt='o')
	if(loop == 4):
		plt.errorbar(daCunha_SM_z4, daCunha_DM_z4, yerr = (daCunha_DMdownerr_z4, daCunha_DMuperr_z4), color='g',label='daCunha2014',fmt='o')
	if(loop == 5):
		plt.errorbar(daCunha_SM_z5, daCunha_DM_z5, yerr = (daCunha_DMdownerr_z5, daCunha_DMuperr_z5), color='g',label='daCunha2014',fmt='o')
	if(loop == 6):
		plt.errorbar(daCunha_SM_z6, daCunha_DM_z6, yerr = (daCunha_DMdownerr_z6, daCunha_DMuperr_z6), color='g',label='daCunha2014',fmt='o')
	if( (loop == 6) or (loop == 7) ):
		plt.errorbar(Mancini_SM, Mancini_DM, yerr = Mancini_DMerr, xerr = Mancini_SMerr, color='g',label='Mancini2015',fmt='o')
	plt.errorbar(SM_bins[4][loop],Dust_bins[4][loop],yerr=(Dust_std_err[0][loop]),color='k',label='tacc=0.2Myr tdes=tdes',linewidth=2)

	#ALL_NO_DEST
	condition = np.logical_and(sSFR_6 > 0.0345E-9, np.logical_and(Dust_Mass_6 > 0, Type_6 == 0))
	log_Stellar_Mass_6 = np.log10(Stellar_Mass_6[condition==1])
	log_New_Dust_Mass_6 = np.log10(Dust_Mass_6[condition==1])
	SM_bins[5][loop],Dust_bins[5][loop],Dust_std_dev[5][loop],Dust_std_err[5][loop],count[5][loop] = fit_scatter(log_Stellar_Mass_6, log_New_Dust_Mass_6, ret_n=True, ret_sterr=True, nbins=10)
	plt.subplot(3,2,5)
	plt.xlim([8,12])
	plt.ylim([0,10.2])
	plt.hexbin(log_Stellar_Mass_6,log_New_Dust_Mass_6,gridsize=500,mincnt=1, label='L-Galaxies 2Dhist')
	plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
	plt.ylabel(r'log$_{10}$(M$_{\rm{d}}$/M$_{\odot}$)', fontsize=18)
	plt.tick_params(axis='both', which='major', labelsize=12,width=2)
	plt.tick_params(axis='both', which='minor', labelsize=12,width=2)
	#plt.text(7,1,"N = "+str(sum(count))+"\nz = "+str(loop))
	plt.text(9.6,1,"z = "+str(loop)+"\nNo Destruction",fontsize=16)
	if(loop == 0):
		plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='g',label='RemyRuyer2014',fmt='o')
		plt.errorbar(np.log10(Bourne_MEDMSTAR), np.log10(Bourne_MEDMDUST), yerr = np.log10(Bourne_MEDMDUSTERR/Bourne_MEDMDUST) , color='orange',label='Bourne2012',fmt='o')
		plt.errorbar(Ciesla_SM, Ciesla_DM, yerr = Ciesla_DMerr , color='r',label='Ciesla2014',fmt='o')
		plt.errorbar(Santini_SM_z0, Santini_DM_z0, yerr = (Santini_DMdownerr_z0, Santini_DMuperr_z0), color='b',label='Santini2014',fmt='o')
	if(loop == 1):
		plt.errorbar(Santini_SM_z1, Santini_DM_z1, yerr = (Santini_DMdownerr_z1, Santini_DMuperr_z1), color='g',label='Santini2014',fmt='o')
	if(loop == 2):
		plt.errorbar(Santini_SM_z2, Santini_DM_z2, yerr = (Santini_DMdownerr_z2, Santini_DMuperr_z2), color='g',label='Santini2014',fmt='o')
		plt.errorbar(daCunha_SM_z2, daCunha_DM_z2, yerr = (daCunha_DMdownerr_z2, daCunha_DMuperr_z2), color='b',label='daCunha2014',fmt='o')
	if(loop == 3):
		plt.errorbar(daCunha_SM_z3, daCunha_DM_z3, yerr = (daCunha_DMdownerr_z3, daCunha_DMuperr_z3), color='g',label='daCunha2014',fmt='o')
	if(loop == 4):
		plt.errorbar(daCunha_SM_z4, daCunha_DM_z4, yerr = (daCunha_DMdownerr_z4, daCunha_DMuperr_z4), color='g',label='daCunha2014',fmt='o')
	if(loop == 5):
		plt.errorbar(daCunha_SM_z5, daCunha_DM_z5, yerr = (daCunha_DMdownerr_z5, daCunha_DMuperr_z5), color='g',label='daCunha2014',fmt='o')
	if(loop == 6):
		plt.errorbar(daCunha_SM_z6, daCunha_DM_z6, yerr = (daCunha_DMdownerr_z6, daCunha_DMuperr_z6), color='g',label='daCunha2014',fmt='o')
	if( (loop == 6) or (loop == 7) ):
		plt.errorbar(Mancini_SM, Mancini_DM, yerr = Mancini_DMerr, xerr = Mancini_SMerr, color='g',label='Mancini2015',fmt='o')
	plt.errorbar(SM_bins[5][loop],Dust_bins[5][loop],yerr=(Dust_std_err[0][loop]),color='k',label='tacc=0.2Myr tdes=tdes',linewidth=2)
	#ALL
	condition = np.logical_and(sSFR_1 > 0.0345E-9, np.logical_and(Dust_Mass_1 > 0, Type_1 == 0))
	log_Stellar_Mass_1 = np.log10(Stellar_Mass_1[condition==1])
	log_New_Dust_Mass_1 = np.log10(Dust_Mass_1[condition==1])
	SM_bins[0][loop],Dust_bins[0][loop],Dust_std_dev[0][loop],Dust_std_err[0][loop],count[0][loop] = fit_scatter(log_Stellar_Mass_1, log_New_Dust_Mass_1, ret_n=True, ret_sterr=True, nbins=10)
	plt.subplot(3,2,6)
	plt.xlim([8,12])
	plt.ylim([0,10.2])
	plt.hexbin(log_Stellar_Mass_1,log_New_Dust_Mass_1,gridsize=500,mincnt=1, label='L-Galaxies 2Dhist')
	plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
	plt.ylabel(r'log$_{10}$(M$_{\rm{d}}$/M$_{\odot}$)', fontsize=18)
	plt.tick_params(axis='both', which='major', labelsize=12,width=2)
	plt.tick_params(axis='both', which='minor', labelsize=12,width=2)
	#plt.text(7,1,"N = "+str(sum(count))+"\nz = "+str(loop))
	plt.text(9.6,1,"z = "+str(loop)+"\nEverything",fontsize=16)
	if(loop == 0):
		plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='g',label='RemyRuyer2014',fmt='o')
		plt.errorbar(np.log10(Bourne_MEDMSTAR), np.log10(Bourne_MEDMDUST), yerr = np.log10(Bourne_MEDMDUSTERR/Bourne_MEDMDUST) , color='orange',label='Bourne2012',fmt='o')
		plt.errorbar(Ciesla_SM, Ciesla_DM, yerr = Ciesla_DMerr , color='r',label='Ciesla2014',fmt='o')
		plt.errorbar(Santini_SM_z0, Santini_DM_z0, yerr = (Santini_DMdownerr_z0, Santini_DMuperr_z0), color='b',label='Santini2014',fmt='o')
	if(loop == 1):
		plt.errorbar(Santini_SM_z1, Santini_DM_z1, yerr = (Santini_DMdownerr_z1, Santini_DMuperr_z1), color='g',label='Santini2014',fmt='o')
	if(loop == 2):
		plt.errorbar(Santini_SM_z2, Santini_DM_z2, yerr = (Santini_DMdownerr_z2, Santini_DMuperr_z2), color='g',label='Santini2014',fmt='o')
		plt.errorbar(daCunha_SM_z2, daCunha_DM_z2, yerr = (daCunha_DMdownerr_z2, daCunha_DMuperr_z2), color='b',label='daCunha2014',fmt='o')
	if(loop == 3):
		plt.errorbar(daCunha_SM_z3, daCunha_DM_z3, yerr = (daCunha_DMdownerr_z3, daCunha_DMuperr_z3), color='g',label='daCunha2014',fmt='o')
	if(loop == 4):
		plt.errorbar(daCunha_SM_z4, daCunha_DM_z4, yerr = (daCunha_DMdownerr_z4, daCunha_DMuperr_z4), color='g',label='daCunha2014',fmt='o')
	if(loop == 5):
		plt.errorbar(daCunha_SM_z5, daCunha_DM_z5, yerr = (daCunha_DMdownerr_z5, daCunha_DMuperr_z5), color='g',label='daCunha2014',fmt='o')
	if(loop == 6):
		plt.errorbar(daCunha_SM_z6, daCunha_DM_z6, yerr = (daCunha_DMdownerr_z6, daCunha_DMuperr_z6), color='g',label='daCunha2014',fmt='o')
	if( (loop == 6) or (loop == 7) ):
		plt.errorbar(Mancini_SM, Mancini_DM, yerr = Mancini_DMerr, xerr = Mancini_SMerr, color='g',label='Mancini2015',fmt='o')
	plt.errorbar(SM_bins[0][loop],Dust_bins[0][loop],yerr=(Dust_std_err[0][loop]),color='k',label='tacc=0.2Myr tdes=tdes',linewidth=2)


	axes = fig.get_axes()
	for ax in axes:
		[i.set_linewidth(2.1) for i in ax.spines.itervalues()]
	
	
	pylab.savefig('./graphs/SM_DM_vary_mechanisms_z'+str(loop)+'.png', bbox_inches=0)
	plt.close()
	
	
	
	
	#---------------------SM vs. DM for varied tdes rates ALL REDSHIFTS
fig = plt.figure(figsize=(9,9))
for loop in range(0,9):
	plt.subplot(3,3,loop+1)
	plt.xlim([8,12])
	plt.ylim([0,10.2])
	plt.tick_params(axis='both', which='major', labelsize=12,width=2)
	plt.tick_params(axis='both', which='minor', labelsize=12,width=2)
	#plt.text(7,1,"N = "+str(sum(count))+"\nz = "+str(loop))
	if(loop !=8):
		plt.text(10,1,"z = "+str(loop),fontsize=16)
	
	if(loop ==8):
		plt.axis('off')
		plt.errorbar([3,4],[3,4],yerr=[3,4],color='b',label='All',linewidth=2)
		plt.errorbar([3,4],[3,4],yerr=[3,4],color='r',label='AGB',linewidth=2)
		plt.errorbar([3,4],[3,4],yerr=[3,4],color='r',linestyle='--',label='SNII',linewidth=2)
		plt.errorbar([3,4],[3,4],yerr=[3,4],color='g',label='SNIA',linewidth=2)	
		plt.errorbar([3,4],[3,4],yerr=[3,4],color='g',linestyle='--',label='MC Grown',linewidth=2)
		plt.errorbar([3,4],[3,4],yerr=[3,4],color='orange',linestyle='--',label='No Destruction',linewidth=2)
		plt.errorbar([3,4],[3,4],yerr=[3,4],color='k',label='Observations',fmt='o')
		plt.legend(loc='center')
	if(loop == 0):
		plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='k',label='RemyRuyer2014',fmt='o')
		plt.errorbar(np.log10(Bourne_MEDMSTAR), np.log10(Bourne_MEDMDUST), yerr = np.log10(Bourne_MEDMDUSTERR/Bourne_MEDMDUST) , color='k',label='Bourne2012',fmt='o')
		plt.errorbar(Ciesla_SM, Ciesla_DM, yerr = Ciesla_DMerr , color='k',label='Ciesla2014',fmt='o')
		plt.errorbar(Santini_SM_z0, Santini_DM_z0, yerr = (Santini_DMdownerr_z0, Santini_DMuperr_z0), color='k',label='Santini2014',fmt='o')
	if(loop == 1):
		plt.errorbar(Santini_SM_z1, Santini_DM_z1, yerr = (Santini_DMdownerr_z1, Santini_DMuperr_z1), color='k',label='Santini2014',fmt='o')
	if(loop == 2):
		plt.errorbar(Santini_SM_z2, Santini_DM_z2, yerr = (Santini_DMdownerr_z2, Santini_DMuperr_z2), color='k',label='Santini2014',fmt='o')
		plt.errorbar(daCunha_SM_z2, daCunha_DM_z2, yerr = (daCunha_DMdownerr_z2, daCunha_DMuperr_z2), color='k',label='daCunha2014',fmt='o')
	if(loop == 3):
		plt.errorbar(daCunha_SM_z3, daCunha_DM_z3, yerr = (daCunha_DMdownerr_z3, daCunha_DMuperr_z3), color='k',label='daCunha2014',fmt='o')
		plt.ylabel(r'log$_{10}$(M$_{\rm{d}}$/M$_{\odot}$)', fontsize=18)
	if(loop == 4):
		plt.errorbar(daCunha_SM_z4, daCunha_DM_z4, yerr = (daCunha_DMdownerr_z4, daCunha_DMuperr_z4), color='k',label='daCunha2014',fmt='o')
	if(loop == 5):
		plt.errorbar(daCunha_SM_z5, daCunha_DM_z5, yerr = (daCunha_DMdownerr_z5, daCunha_DMuperr_z5), color='k',label='daCunha2014',fmt='o')
	if(loop == 6):
		plt.errorbar(daCunha_SM_z6, daCunha_DM_z6, yerr = (daCunha_DMdownerr_z6, daCunha_DMuperr_z6), color='k',label='daCunha2014',fmt='o')
	if( (loop == 6) or (loop == 7) ):
		plt.errorbar(Mancini_SM, Mancini_DM, yerr = Mancini_DMerr, xerr = Mancini_SMerr, color='k',label='Mancini2015',fmt='o')
	if(loop == 7):
		plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
	if(loop !=8):
		plt.errorbar(SM_bins[0][loop],Dust_bins[0][loop],yerr=(Dust_std_err[0][loop]),color='b',label='L-Galaxies Mean',linewidth=2)
		plt.errorbar(SM_bins[1][loop],Dust_bins[1][loop],yerr=(Dust_std_err[1][loop]),color='r',label='L-Galaxies Mean',linewidth=2)
		plt.errorbar(SM_bins[2][loop],Dust_bins[2][loop],yerr=(Dust_std_err[2][loop]),color='r',linestyle='--',label='L-Galaxies Mean',linewidth=2)
		plt.errorbar(SM_bins[3][loop],Dust_bins[3][loop],yerr=(Dust_std_err[3][loop]),color='g',label='L-Galaxies Mean',linewidth=2)	
		plt.errorbar(SM_bins[4][loop],Dust_bins[4][loop],yerr=(Dust_std_err[4][loop]),color='g',linestyle='--',label='L-Galaxies Mean',linewidth=2)
		plt.errorbar(SM_bins[5][loop],Dust_bins[5][loop],yerr=(Dust_std_err[5][loop]),color='orange',linestyle='--',label='L-Galaxies Mean',linewidth=2)

# plt.legend(loc='lower right')

axes = fig.get_axes()
for ax in axes:
	[i.set_linewidth(2.1) for i in ax.spines.itervalues()]


pylab.savefig('./graphs/SM_DM_vary_mechanisms_ALLz.png', bbox_inches=0)
plt.close()
