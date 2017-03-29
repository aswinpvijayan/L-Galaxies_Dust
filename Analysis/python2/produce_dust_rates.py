import cPickle
import numpy as np
import matplotlib.pyplot as plt
import pylab
import sys


DustRate_z = np.zeros(10)
DustRate_AGB_z  = np.zeros(10)
DustRate_SNII_z = np.zeros(10)
DustRate_SNIA_z = np.zeros(10)
DustRate_GROW_z = np.zeros(10)
DustRate_DEST_z = np.zeros(10)
DustRate_ALL_z = np.zeros(10)

for loop in range(0,10):

	fin = open('../data/'+str(sys.argv[1])+'/lgal_z'+str(loop)+'.pkl','rb')
	gals=cPickle.load(fin)
	fin.close()

	Type = np.zeros(len(gals['Type']))
	Stellar_Mass = np.zeros(len(gals['Type']))
	Stellar_Mass = np.log10(gals['StellarMass']*1.0E10/0.673)

	DustRate_AGB  = np.zeros(len(gals['Type']))
	DustRate_SNII = np.zeros(len(gals['Type']))
	DustRate_SNIA = np.zeros(len(gals['Type']))
	DustRate_GROW = np.zeros(len(gals['Type']))
	DustRate_DEST = np.zeros(len(gals['Type']))
	DustRate_ALL = np.zeros(len(gals['Type']))

	for i in range(0,len(gals['Type'])):
		DustRate_AGB[i]  = gals['DustRatesISM'][i][0]
		DustRate_SNII[i] = gals['DustRatesISM'][i][1]
		DustRate_SNIA[i] = gals['DustRatesISM'][i][2]
		DustRate_GROW[i] = gals['DustRatesISM'][i][3]
		DustRate_DEST[i] = gals['DustRatesISM'][i][4]
		DustRate_ALL[i] = gals['DustRatesISM'][i][0]+gals['DustRatesISM'][i][1]+gals['DustRatesISM'][i][2]+gals['DustRatesISM'][i][3]

	DustRate_z[loop] = loop

	condition = np.logical_and(Stellar_Mass > (6.0),DustRate_AGB>0.0)
	DustRate_AGB_z[loop]= np.median(np.log10(DustRate_AGB[condition==1]))

	condition = np.logical_and(Stellar_Mass > (6.0),DustRate_SNII>0.0)
	DustRate_SNII_z[loop]= np.median(np.log10(DustRate_SNII[condition==1]))

	condition = np.logical_and(Stellar_Mass > (6.0),DustRate_SNIA>0.0)
	DustRate_SNIA_z[loop]= np.median(np.log10(DustRate_SNIA[condition==1]))

	condition = np.logical_and(Stellar_Mass > (6.0),DustRate_GROW>0.0)
	DustRate_GROW_z[loop]= np.median(np.log10(DustRate_GROW[condition==1]))

	condition = np.logical_and(Stellar_Mass > (6.0),DustRate_DEST>0.0)
	DustRate_DEST_z[loop]= np.median(np.log10(DustRate_DEST[condition==1]))

	condition = np.logical_and(Stellar_Mass > (6.0),DustRate_ALL>0.0)
	DustRate_ALL_z[loop]= np.median(np.log10(DustRate_ALL[condition==1]))

plt.xlabel(r'redshift', fontsize=14,labelpad=10)
plt.ylabel(r'log$_{10}$(Cosmic dust rate Msol/yr)}$)', fontsize=14,labelpad=0)
plt.xlim(-1,12)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.tick_params(axis='both', which='minor', labelsize=8)

plt.plot(DustRate_z,DustRate_AGB_z , color='r',label='AGB')
plt.plot(DustRate_z,DustRate_SNII_z, color='b',label='SNII')
plt.plot(DustRate_z,DustRate_SNIA_z, color='y',label='SNIA')
plt.plot(DustRate_z,DustRate_GROW_z, color='g',label='GROW')
plt.plot(DustRate_z,DustRate_DEST_z, color='k',label='DEST')
plt.plot(DustRate_z,DustRate_ALL_z,  color='orange',label='ALL')
plt.legend(loc='lower right')
pylab.savefig('./graphs/dust_redshift_ALLMASS.png', bbox_inches=0)
plt.close()

for mass in [9.0,10.0,11.0,12.0]:
	DustRate_z = np.zeros(10)
	DustRate_AGB_z  = np.zeros(10)
	DustRate_SNII_z = np.zeros(10)
	DustRate_SNIA_z = np.zeros(10)
	DustRate_GROW_z = np.zeros(10)
	DustRate_DEST_z = np.zeros(10)
	DustRate_ALL_z = np.zeros(10)

	for loop in range(0,10):
	
		fin = open('../data/'+str(sys.argv[1])+'/lgal_z'+str(loop)+'.pkl','rb')
		gals=cPickle.load(fin)
		fin.close()

		Type = np.zeros(len(gals['Type']))
		Stellar_Mass = np.zeros(len(gals['Type']))
		Stellar_Mass = np.log10(gals['StellarMass']*1.0E10/0.673)

		DustRate_AGB  = np.zeros(len(gals['Type']))
		DustRate_SNII = np.zeros(len(gals['Type']))
		DustRate_SNIA = np.zeros(len(gals['Type']))
		DustRate_GROW = np.zeros(len(gals['Type']))
		DustRate_DEST = np.zeros(len(gals['Type']))
		DustRate_ALL = np.zeros(len(gals['Type']))

		for i in range(0,len(gals['Type'])):
			DustRate_AGB[i]  = gals['DustRatesISM'][i][0]
			DustRate_SNII[i] = gals['DustRatesISM'][i][1]
			DustRate_SNIA[i] = gals['DustRatesISM'][i][2]
			DustRate_GROW[i] = gals['DustRatesISM'][i][3]
			DustRate_DEST[i] = gals['DustRatesISM'][i][4]
			DustRate_ALL[i]  = gals['DustRatesISM'][i][0] + gals['DustRatesISM'][i][1]+gals['DustRatesISM'][i][2] + gals['DustRatesISM'][i][3]
	
		DustRate_z[loop] = loop

		condition = np.logical_and(np.logical_and(Stellar_Mass > (mass-0.5),Stellar_Mass < (mass+0.5)),DustRate_AGB>0.0)
		DustRate_AGB_z[loop]= np.median(np.log10(DustRate_AGB[condition==1]))

		condition = np.logical_and(np.logical_and(Stellar_Mass > (mass-0.5),Stellar_Mass < (mass+0.5)),DustRate_SNII>0.0)
		DustRate_SNII_z[loop]= np.median(np.log10(DustRate_SNII[condition==1]))

		condition = np.logical_and(np.logical_and(Stellar_Mass > (mass-0.5),Stellar_Mass < (mass+0.5)),DustRate_SNIA>0.0)
		DustRate_SNIA_z[loop]= np.median(np.log10(DustRate_SNIA[condition==1]))

		condition = np.logical_and(np.logical_and(Stellar_Mass > (mass-0.5),Stellar_Mass < (mass+0.5)),DustRate_GROW>0.0)
		DustRate_GROW_z[loop]= np.median(np.log10(DustRate_GROW[condition==1]))

		condition = np.logical_and(np.logical_and(Stellar_Mass > (mass-0.5),Stellar_Mass < (mass+0.5)),DustRate_DEST>0.0)
		DustRate_DEST_z[loop]= np.median(np.log10(DustRate_DEST[condition==1]))

		condition = np.logical_and(np.logical_and(Stellar_Mass > (mass-0.5),Stellar_Mass < (mass+0.5)),DustRate_ALL>0.0)
		DustRate_ALL_z[loop]= np.median(np.log10(DustRate_ALL[condition==1]))



	plt.xlabel(r'redshift', fontsize=14,labelpad=10)
	plt.ylabel(r'log$_{10}$(Cosmic dust rate Msol/yr)}$)', fontsize=14,labelpad=0)
	plt.xlim(-1,12)
	plt.tick_params(axis='both', which='major', labelsize=10)
	plt.tick_params(axis='both', which='minor', labelsize=8)
	
	plt.plot(DustRate_z,DustRate_AGB_z , color='r',label='AGB')
	plt.plot(DustRate_z,DustRate_SNII_z, color='b',label='SNII')
	plt.plot(DustRate_z,DustRate_SNIA_z, color='y',label='SNIA')
	plt.plot(DustRate_z,DustRate_GROW_z, color='g',label='GROW')
	plt.plot(DustRate_z,DustRate_DEST_z, color='k',label='DEST')
	plt.plot(DustRate_z,DustRate_ALL_z, color='orange',label='ALL')
	
# 	plt.errorbar(DustRate_z,avg_AGB_dust ,std_AGB_dust , color='r',label='AGB')
# 	plt.errorbar(DustRate_z,avg_SNII_dust,std_SNII_dust, color='b',label='SNII')
# 	plt.errorbar(DustRate_z,avg_GROW_dust,std_GROW_dust, color='g',label='Grown')
# 	plt.errorbar(DustRate_z,avg_SNIa_dust,std_SNIa_dust, color='y',label='SNIA')
# 	plt.errorbar(DustRate_z,avg_ALL_dust,std_ALL_dust  , color='k',label='ALL')
	plt.legend(loc='lower right')
	pylab.savefig('./graphs/dust_redshift'+str(mass)+'.png', bbox_inches=0)
	plt.close()

