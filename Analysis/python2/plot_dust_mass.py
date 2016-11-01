# Script to read in pickled file of galaxy properties.
# Requires pickleFile to hold name of desired input file.
# Returns galaxy properties in gals.

import cPickle
import numpy as np
import matplotlib.pyplot as plt
import pylab

for loop in range(0,10):
	fin = open('../data/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()


	AGB_Dust_Mass = np.zeros(len(gals['Type']))
	SNII_Dust_Mass = np.zeros(len(gals['Type']))
	SNIa_Dust_Mass = np.zeros(len(gals['Type']))
	Growth_Dust_Mass = np.zeros(len(gals['Type']))
	Stellar_Mass = np.zeros(len(gals['Type']))

	for i in range(0,len(gals['Type'])):
		AGB_Dust_Mass[i] = (gals['DustMassISM'][i][0]+gals['DustMassISM'][i][1]+gals['DustMassISM'][i][2]+gals['DustMassISM'][i][3])
		SNII_Dust_Mass[i] = (gals['DustMassISM'][i][4]+gals['DustMassISM'][i][5]+gals['DustMassISM'][i][6]+gals['DustMassISM'][i][7])
		SNIa_Dust_Mass[i] = (gals['DustMassISM'][i][8]+gals['DustMassISM'][i][9]+gals['DustMassISM'][i][10]+gals['DustMassISM'][i][11])
		Growth_Dust_Mass[i] = (gals['DustMassISM'][i][12]+gals['DustMassISM'][i][13]+gals['DustMassISM'][i][14]+gals['DustMassISM'][i][15])
		Stellar_Mass[i] = gals['StellarMass'][i]*1.0E10/0.673

	log_AGB = np.log10(AGB_Dust_Mass)
	log_SM = np.log10(Stellar_Mass)


	# Mean values
	avg=np.zeros(5)
	stdev=np.zeros(5)

	j=0
	xxx=[7,8,9,10,11]
	for i in range(7,12):
		low=i
		high=low+1
	#	avg[j]=np.median(AGB_Dust_Mass[(np.log10(Stellar_Mass)>=low)&(np.log10(Stellar_Mass)<high)])
	#	stdev[j] = np.std(AGB_Dust_Mass[(np.log10(Stellar_Mass)>=low)&(np.log10(Stellar_Mass)<high)])
		avg[j]=np.median(log_AGB[(log_SM>=low)&(log_SM<high)])
		stdev[j] =np.std(log_AGB[(log_SM>=low)&(log_SM<high)])

		print i, avg[j],stdev[j]
		j=j+1
	


	# Plot 2D histogram using pcolor
	#plt.figure(figsize=(5,5))
	#plt.ylim([-10,-7])
	plt.xlim([7,12])
	#plt.pcolormesh(xedges,yedges,Hmasked)
	plt.hexbin(log_SM,log_AGB,gridsize=500,mincnt=1)
	#plt.plot(xxx,np.log10(avg),color='r')
	plt.errorbar(xxx,(avg),yerr=(stdev),color='r')

	plt.xlabel(r'log$_{10}$(Mstellar/M$_{\odot}$)', fontsize=14,labelpad=10)
	plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
	plt.tick_params(axis='both', which='major', labelsize=10)
	plt.tick_params(axis='both', which='minor', labelsize=8)
	#cbar = plt.colorbar()
	#cbar.ax.set_ylabel('Counts',fontsize=14)
	#cbar.ax.tick_params(labelsize=6) 


	#plt.plot(MASS,METALLICITY,'ro')

	#plt.show()
	pylab.savefig('stellar_dust_z'+str(loop)+'.png', bbox_inches=0)
	plt.close()









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

















