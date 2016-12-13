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

obs_DM, obs_DM_err, obs_SM = np.loadtxt('../observations/z0.txt',unpack=True)
obs_DM_act_err = obs_DM*(obs_DM_err/100.0)

print obs_DM

for loop in range(0,10):
#for loop in range(0,1):

	fin = open('../data/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()


	AGB_Dust_Mass = np.zeros(len(gals['Type']))
	SNII_Dust_Mass = np.zeros(len(gals['Type']))
	SNIa_Dust_Mass = np.zeros(len(gals['Type']))
	Growth_Dust_Mass = np.zeros(len(gals['Type']))
	Stellar_Mass = np.zeros(len(gals['Type']))
	All_Dust_Mass = np.zeros(len(gals['Type']))
	
	AGB_Dust_SiC= np.zeros(len(gals['Type']))
	AGB_Dust_Sil= np.zeros(len(gals['Type']))
	AGB_Dust_Fe= np.zeros(len(gals['Type']))
	AGB_Dust_Cb= np.zeros(len(gals['Type']))


	for i in range(0,len(gals['Type'])):
		AGB_Dust_Mass[i] = (gals['DustMassISM'][i][0]+gals['DustMassISM'][i][1]+gals['DustMassISM'][i][2]+gals['DustMassISM'][i][3])
		SNII_Dust_Mass[i] = (gals['DustMassISM'][i][4]+gals['DustMassISM'][i][5]+gals['DustMassISM'][i][6]+gals['DustMassISM'][i][7])
		SNIa_Dust_Mass[i] = (gals['DustMassISM'][i][8]+gals['DustMassISM'][i][9]+gals['DustMassISM'][i][10]+gals['DustMassISM'][i][11])
		Growth_Dust_Mass[i] = (gals['DustMassISM'][i][12]+gals['DustMassISM'][i][13]+gals['DustMassISM'][i][14]+gals['DustMassISM'][i][15])
		Stellar_Mass[i] = gals['StellarMass'][i]*1.0E10/0.673
		
		#print AGB_Dust_Mass[i], SNII_Dust_Mass[i]
		if AGB_Dust_Mass[i] < 1.00 :
			print AGB_Dust_Mass[i]
		
		AGB_Dust_SiC[i] = gals['DustMassISM'][i][0]
		AGB_Dust_Sil[i] = gals['DustMassISM'][i][1]
		AGB_Dust_Fe[i]  = gals['DustMassISM'][i][2]
		AGB_Dust_Cb[i]  = gals['DustMassISM'][i][3]
		
		for j in range(0,16):
			All_Dust_Mass[i] += gals['DustMassISM'][i][j]
			
		

		
		
		
		
# 	log_AGB = np.log10(AGB_Dust_Mass)
	log_SM = np.log10(Stellar_Mass)
	
# 	
# 	
# 	log_AGB_SiC = np.log10(AGB_Dust_SiC)
# 	log_AGB_Sil = np.log10(AGB_Dust_Sil)
# 	log_AGB_Cb  = np.log10(AGB_Dust_Cb)
# 	log_AGB_Fe  = np.log10(AGB_Dust_Fe)
# 
# 	log_SNII = np.log10(SNII_Dust_Mass)
# 
# 	log_GROW = np.log10(Growth_Dust_Mass)
# 
# 
# 	
# 	log_ALL = np.log10(All_Dust_Mass)



	#---------------------AGB
	# Mean values
	avg_AGB=np.zeros(5)
	avg_AGB_SiC=np.zeros(5)
	avg_AGB_Sil=np.zeros(5)
	avg_AGB_Fe=np.zeros(5)
	avg_AGB_Cb=np.zeros(5)
	
	stdev_AGB=np.zeros(5)
	stdev_AGB_SiC=np.zeros(5)
	stdev_AGB_Sil=np.zeros(5)
	stdev_AGB_Fe=np.zeros(5)
	stdev_AGB_Cb=np.zeros(5)
	
	

	j=0
	xxx=[7.5,8.5,9.5,10.5,11.5]
	for i in range(7,12):
		low=i
		high=low+1
		avg_AGB[j]=np.median(log_AGB[(log_SM>=low)&(log_SM<high)])
		stdev_AGB[j] =np.std(log_AGB[(log_SM>=low)&(log_SM<high)])
		avg_AGB_SiC[j]=np.median(log_AGB_SiC[(log_SM>=low)&(log_SM<high)])
		stdev_AGB_SiC[j] =np.std(log_AGB_SiC[(log_SM>=low)&(log_SM<high)])
		avg_AGB_Sil[j]=np.median(log_AGB_Sil[(log_SM>=low)&(log_SM<high)])
		stdev_AGB_Sil[j] =np.std(log_AGB_Sil[(log_SM>=low)&(log_SM<high)])
		avg_AGB_Fe[j]=np.median(log_AGB_Fe[(log_SM>=low)&(log_SM<high)])
		stdev_AGB_Fe[j] =np.std(log_AGB_Fe[(log_SM>=low)&(log_SM<high)])
		avg_AGB_Cb[j]=np.median(log_AGB_Cb[(log_SM>=low)&(log_SM<high)])
		stdev_AGB_Cb[j] =np.std(log_AGB_Cb[(log_SM>=low)&(log_SM<high)])
		#print i, avg[j],stdev[j]
		j=j+1
	
	xx,yy,zz = fit_scatter(Stellar_Mass, AGB_Dust_Mass)
	print xx
	print yy
	print zz
	
	
	plt.xlim([6,12])
	plt.hexbin(log_SM,log_AGB,gridsize=500,mincnt=1, label='Dust AGB')
	plt.errorbar(xxx,(avg_AGB),yerr=(stdev_AGB),color='r',label='Dust AGB')
	plt.errorbar(xx,yy,yerr=zz,color='g',label='Dave')

# 	plt.errorbar(xxx,(avg_AGB_SiC),yerr=(stdev_AGB_SiC),color='g')
# 	plt.errorbar(xxx,(avg_AGB_Sil),yerr=(stdev_AGB_Sil),color='g')
# 	plt.errorbar(xxx,(avg_AGB_Fe),yerr=(stdev_AGB_Fe),color='g')
# 	plt.errorbar(xxx,(avg_AGB_Cb),yerr=(stdev_AGB_Cb),color='g')

	plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
	plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
	plt.tick_params(axis='both', which='major', labelsize=10)
	plt.tick_params(axis='both', which='minor', labelsize=8)
	plt.legend(loc='lower right')
	pylab.savefig('./graphs/stellar_AGBdust_z'+str(loop)+'.png', bbox_inches=0)
	plt.close()
# 
# 	#---------------------SNII
# 	# Mean values
# 	avg_SNII=np.zeros(5)
# 	stdev_SNII=np.zeros(5)
# 
# 	j=0
# 	xxx=[7.5,8.5,9.5,10.5,11.5]
# 	for i in range(7,12):
# 		low=i
# 		high=low+1
# 		avg_SNII[j]=np.median(log_SNII[(log_SM>=low)&(log_SM<high)])
# 		stdev_SNII[j] =np.std(log_SNII[(log_SM>=low)&(log_SM<high)])
# 		#print i, avg[j],stdev[j]
# 		j=j+1
# 	
# 	plt.xlim([6,12])
# 	plt.hexbin(log_SM,log_SNII,gridsize=500,mincnt=1,label='Dust SNII')
# 	plt.errorbar(xxx,(avg_SNII),yerr=(stdev_SNII),color='r',label='Dust SNII')
# 	plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
# 	plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
# 	plt.tick_params(axis='both', which='major', labelsize=10)
# 	plt.tick_params(axis='both', which='minor', labelsize=8)
# 	plt.legend(loc='lower right')
# 	pylab.savefig('./graphs/stellar_SNIIdust_z'+str(loop)+'.png', bbox_inches=0)
# 	plt.close()
# 
# 
# 	#---------------------GROWTH
# 	# Mean values
# 	avg_GROW=np.zeros(5)
# 	stdev_GROW=np.zeros(5)
# 
# 	j=0
# 	xxx=[7.5,8.5,9.5,10.5,11.5]
# 	for i in range(7,12):
# 		low=i
# 		high=low+1
# 		avg_GROW[j]=np.median(log_GROW[(log_SM>=low)&(log_SM<high)])
# 		stdev_GROW[j] =np.std(log_GROW[(log_SM>=low)&(log_SM<high)])
# 		#print i, avg[j],stdev[j]
# 		j=j+1
# 	
# 	plt.xlim([6,12])
# 	plt.hexbin(log_SM,log_GROW,gridsize=500,mincnt=1,label='Dust Growth')
# 	plt.errorbar(xxx,(avg_GROW),yerr=(stdev_GROW),color='r',label='Dust Growth')
# 	plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
# 	plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
# 	plt.tick_params(axis='both', which='major', labelsize=10)
# 	plt.tick_params(axis='both', which='minor', labelsize=8)
# 	plt.legend(loc='lower right')
# 	pylab.savefig('./graphs/stellar_GROWdust_z'+str(loop)+'.png', bbox_inches=0)
# 	plt.close()

# 	---------------------ALL
# 	Mean values
# 	avg_ALL=np.zeros(5)
# 	stdev_ALL=np.zeros(5)
# 
# 	j=0
# 	xxx=[7.5,8.5,9.5,10.5,11.5]
# 	for i in range(7,12):
# 		low=i
# 		high=low+1
# 		avg_ALL[j]=np.mean(log_ALL[(log_SM>=low)&(log_SM<high)])
# 		stdev_ALL[j] =np.std(log_ALL[(log_SM>=low)&(log_SM<high)])
# 		print i, avg[j],stdev[j]
# 		j=j+1
# 	
# 	plt.xlim([6,12])
# 	plt.hexbin(log_SM,log_ALL,gridsize=500,mincnt=1,label='Dust ALL')
# 	plt.errorbar(xxx,(avg_ALL),yerr=(stdev_ALL),color='r',label='Dust ALL')
# 	if(loop == 0):
# 		plt.scatter(obs_SM,np.log10(obs_DM),yerr=np.log10(obs_DM_err),color='g',label='RemyRuyer2014')
# 		plt.errorbar(obs_SM,np.log10(obs_DM),yerr=np.log10(obs_DM)*(obs_DM_err/100.0),color='g',label='RemyRuyer2014',fmt='o')
# 
# 	plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
# 	plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
# 	plt.tick_params(axis='both', which='major', labelsize=10)
# 	plt.tick_params(axis='both', which='minor', labelsize=8)
# 	plt.legend(loc='lower right')
# 	pylab.savefig('./graphs/stellar_ALLdust_z'+str(loop)+'.png', bbox_inches=0)
# 	plt.close()


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

















