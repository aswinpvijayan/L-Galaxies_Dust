# Script to read in pickled file of galaxy properties.
# Requires pickleFile to hold name of desired input file.
# Returns galaxy properties in gals.

import cPickle
import numpy as np
import matplotlib.pyplot as plt
import pylab

fin = open('../data/lgal_z0.pkl','rb')
gals_z0=cPickle.load(fin)
fin.close()
fin = open('../data/lgal_z1.pkl','rb')
gals_z1=cPickle.load(fin)
fin.close()
fin = open('../data/lgal_z2.pkl','rb')
gals_z2=cPickle.load(fin)
fin.close()
fin = open('../data/lgal_z3.pkl','rb')
gals_z3=cPickle.load(fin)
fin.close()
fin = open('../data/lgal_z4.pkl','rb')
gals_z4=cPickle.load(fin)
fin.close()
fin = open('../data/lgal_z5.pkl','rb')
gals_z5=cPickle.load(fin)
fin.close()
fin = open('../data/lgal_z6.pkl','rb')
gals_z6=cPickle.load(fin)
fin.close()
fin = open('../data/lgal_z7.pkl','rb')
gals_z7=cPickle.load(fin)
fin.close()
fin = open('../data/lgal_z8.pkl','rb')
gals_z8=cPickle.load(fin)
fin.close()
fin = open('../data/lgal_z9.pkl','rb')
gals_z9=cPickle.load(fin)
fin.close()


AGB_Dust_Mass_z0 = np.zeros(len(gals_z0['Type']))
SNII_Dust_Mass_z0 = np.zeros(len(gals_z0['Type']))
SNIA_Dust_Mass_z0 = np.zeros(len(gals_z0['Type']))
Growth_Dust_Mass_z0 = np.zeros(len(gals_z0['Type']))
Stellar_Mass_z0 = np.zeros(len(gals_z0['Type']))

AGB_Dust_Mass_z1 = np.zeros(len(gals_z1['Type']))
SNII_Dust_Mass_z1 = np.zeros(len(gals_z1['Type']))
SNIA_Dust_Mass_z1 = np.zeros(len(gals_z1['Type']))
Growth_Dust_Mass_z1 = np.zeros(len(gals_z1['Type']))
Stellar_Mass_z1 = np.zeros(len(gals_z1['Type']))

AGB_Dust_Mass_z2 = np.zeros(len(gals_z2['Type']))
SNII_Dust_Mass_z2 = np.zeros(len(gals_z2['Type']))
SNIA_Dust_Mass_z2 = np.zeros(len(gals_z2['Type']))
Growth_Dust_Mass_z2 = np.zeros(len(gals_z2['Type']))
Stellar_Mass_z2 = np.zeros(len(gals_z2['Type']))

AGB_Dust_Mass_z3 = np.zeros(len(gals_z3['Type']))
SNII_Dust_Mass_z3 = np.zeros(len(gals_z3['Type']))
SNIA_Dust_Mass_z3 = np.zeros(len(gals_z3['Type']))
Growth_Dust_Mass_z3 = np.zeros(len(gals_z3['Type']))
Stellar_Mass_z3 = np.zeros(len(gals_z3['Type']))

AGB_Dust_Mass_z4 = np.zeros(len(gals_z4['Type']))
SNII_Dust_Mass_z4 = np.zeros(len(gals_z4['Type']))
SNIA_Dust_Mass_z4 = np.zeros(len(gals_z4['Type']))
Growth_Dust_Mass_z4 = np.zeros(len(gals_z4['Type']))
Stellar_Mass_z4 = np.zeros(len(gals_z4['Type']))

AGB_Dust_Mass_z5 = np.zeros(len(gals_z5['Type']))
SNII_Dust_Mass_z5 = np.zeros(len(gals_z5['Type']))
SNIA_Dust_Mass_z5 = np.zeros(len(gals_z5['Type']))
Growth_Dust_Mass_z5 = np.zeros(len(gals_z5['Type']))
Stellar_Mass_z5 = np.zeros(len(gals_z5['Type']))

AGB_Dust_Mass_z6 = np.zeros(len(gals_z6['Type']))
SNII_Dust_Mass_z6 = np.zeros(len(gals_z6['Type']))
SNIA_Dust_Mass_z6 = np.zeros(len(gals_z6['Type']))
Growth_Dust_Mass_z6 = np.zeros(len(gals_z6['Type']))
Stellar_Mass_z6 = np.zeros(len(gals_z6['Type']))

AGB_Dust_Mass_z7 = np.zeros(len(gals_z7['Type']))
SNII_Dust_Mass_z7 = np.zeros(len(gals_z7['Type']))
SNIA_Dust_Mass_z7 = np.zeros(len(gals_z7['Type']))
Growth_Dust_Mass_z7 = np.zeros(len(gals_z7['Type']))
Stellar_Mass_z7 = np.zeros(len(gals_z7['Type']))

AGB_Dust_Mass_z8 = np.zeros(len(gals_z8['Type']))
SNII_Dust_Mass_z8 = np.zeros(len(gals_z8['Type']))
SNIA_Dust_Mass_z8 = np.zeros(len(gals_z8['Type']))
Growth_Dust_Mass_z8 = np.zeros(len(gals_z8['Type']))
Stellar_Mass_z8 = np.zeros(len(gals_z8['Type']))

AGB_Dust_Mass_z9 = np.zeros(len(gals_z9['Type']))
SNII_Dust_Mass_z9 = np.zeros(len(gals_z9['Type']))
SNIA_Dust_Mass_z9 = np.zeros(len(gals_z9['Type']))
Growth_Dust_Mass_z9 = np.zeros(len(gals_z9['Type']))
Stellar_Mass_z9 = np.zeros(len(gals_z9['Type']))

for i in range(0,len(gals_z0['Type'])):
	AGB_Dust_Mass_z0[i]  = np.log10(gals_z0['DustMassISM'][i][0]+gals_z0['DustMassISM'][i][1]+gals_z0['DustMassISM'][i][2]+gals_z0['DustMassISM'][i][3])
	SNII_Dust_Mass_z0[i] = np.log10(gals_z0['DustMassISM'][i][4]+gals_z0['DustMassISM'][i][5]+gals_z0['DustMassISM'][i][6]+gals_z0['DustMassISM'][i][7])
	SNIA_Dust_Mass_z0[i] = np.log10(gals_z0['DustMassISM'][i][8]+gals_z0['DustMassISM'][i][9]+gals_z0['DustMassISM'][i][10]+gals_z0['DustMassISM'][i][11])
	Growth_Dust_Mass_z0[i] = np.log10(gals_z0['DustMassISM'][i][12]+gals_z0['DustMassISM'][i][13]+gals_z0['DustMassISM'][i][14]+gals_z0['DustMassISM'][i][15])
	Stellar_Mass_z0[i] = np.log10(gals_z0['StellarMass'][i]*1.0E10/0.673)

for i in range(0,len(gals_z1['Type'])):
	AGB_Dust_Mass_z1[i] = np.log10(gals_z1['DustMassISM'][i][0]+gals_z1['DustMassISM'][i][1]+gals_z1['DustMassISM'][i][2]+gals_z1['DustMassISM'][i][3])
	SNII_Dust_Mass_z1[i] = np.log10(gals_z1['DustMassISM'][i][4]+gals_z1['DustMassISM'][i][5]+gals_z1['DustMassISM'][i][6]+gals_z1['DustMassISM'][i][7])
	SNIA_Dust_Mass_z1[i] = np.log10(gals_z1['DustMassISM'][i][8]+gals_z1['DustMassISM'][i][9]+gals_z1['DustMassISM'][i][10]+gals_z1['DustMassISM'][i][11])
	Growth_Dust_Mass_z1[i] = np.log10(gals_z1['DustMassISM'][i][12]+gals_z1['DustMassISM'][i][13]+gals_z1['DustMassISM'][i][14]+gals_z1['DustMassISM'][i][15])
	Stellar_Mass_z1[i] = np.log10(gals_z1['StellarMass'][i]*1.0E10/0.673)

for i in range(0,len(gals_z2['Type'])):
	AGB_Dust_Mass_z2[i] = np.log10(gals_z2['DustMassISM'][i][0]+gals_z2['DustMassISM'][i][1]+gals_z2['DustMassISM'][i][2]+gals_z2['DustMassISM'][i][3])
	SNII_Dust_Mass_z2[i] =np.log10 (gals_z2['DustMassISM'][i][4]+gals_z2['DustMassISM'][i][5]+gals_z2['DustMassISM'][i][6]+gals_z2['DustMassISM'][i][7])
	SNIA_Dust_Mass_z2[i] =np.log10 (gals_z2['DustMassISM'][i][8]+gals_z2['DustMassISM'][i][9]+gals_z2['DustMassISM'][i][10]+gals_z2['DustMassISM'][i][11])
	Growth_Dust_Mass_z2[i] = np.log10(gals_z2['DustMassISM'][i][12]+gals_z2['DustMassISM'][i][13]+gals_z2['DustMassISM'][i][14]+gals_z2['DustMassISM'][i][15])
	Stellar_Mass_z2[i] = np.log10(gals_z2['StellarMass'][i]*1.0E10/0.673)

for i in range(0,len(gals_z3['Type'])):
	AGB_Dust_Mass_z3[i] = np.log10(gals_z3['DustMassISM'][i][0]+gals_z3['DustMassISM'][i][1]+gals_z3['DustMassISM'][i][2]+gals_z3['DustMassISM'][i][3])
	SNII_Dust_Mass_z3[i] =np.log10 (gals_z3['DustMassISM'][i][4]+gals_z3['DustMassISM'][i][5]+gals_z3['DustMassISM'][i][6]+gals_z3['DustMassISM'][i][7])
	SNIA_Dust_Mass_z3[i] = np.log10(gals_z3['DustMassISM'][i][8]+gals_z3['DustMassISM'][i][9]+gals_z3['DustMassISM'][i][10]+gals_z3['DustMassISM'][i][11])
	Growth_Dust_Mass_z3[i] =np.log10 (gals_z3['DustMassISM'][i][12]+gals_z3['DustMassISM'][i][13]+gals_z3['DustMassISM'][i][14]+gals_z3['DustMassISM'][i][15])
	Stellar_Mass_z3[i] = np.log10(gals_z3['StellarMass'][i]*1.0E10/0.673)

for i in range(0,len(gals_z4['Type'])):
	AGB_Dust_Mass_z4[i] =np.log10 (gals_z4['DustMassISM'][i][0]+gals_z4['DustMassISM'][i][1]+gals_z4['DustMassISM'][i][2]+gals_z4['DustMassISM'][i][3])
	SNII_Dust_Mass_z4[i] = np.log10(gals_z4['DustMassISM'][i][4]+gals_z4['DustMassISM'][i][5]+gals_z4['DustMassISM'][i][6]+gals_z4['DustMassISM'][i][7])
	SNIA_Dust_Mass_z4[i] = np.log10(gals_z4['DustMassISM'][i][8]+gals_z4['DustMassISM'][i][9]+gals_z4['DustMassISM'][i][10]+gals_z4['DustMassISM'][i][11])
	Growth_Dust_Mass_z4[i] = np.log10(gals_z4['DustMassISM'][i][12]+gals_z4['DustMassISM'][i][13]+gals_z4['DustMassISM'][i][14]+gals_z4['DustMassISM'][i][15])
	Stellar_Mass_z4[i] = np.log10(gals_z4['StellarMass'][i]*1.0E10/0.673)

for i in range(0,len(gals_z5['Type'])):
	AGB_Dust_Mass_z5[i] =np.log10 (gals_z5['DustMassISM'][i][0]+gals_z5['DustMassISM'][i][1]+gals_z5['DustMassISM'][i][2]+gals_z5['DustMassISM'][i][3])
	SNII_Dust_Mass_z5[i] = np.log10(gals_z5['DustMassISM'][i][4]+gals_z5['DustMassISM'][i][5]+gals_z5['DustMassISM'][i][6]+gals_z5['DustMassISM'][i][7])
	SNIA_Dust_Mass_z5[i] = np.log10(gals_z5['DustMassISM'][i][8]+gals_z5['DustMassISM'][i][9]+gals_z5['DustMassISM'][i][10]+gals_z5['DustMassISM'][i][11])
	Growth_Dust_Mass_z5[i] =np.log10 (gals_z5['DustMassISM'][i][12]+gals_z5['DustMassISM'][i][13]+gals_z5['DustMassISM'][i][14]+gals_z5['DustMassISM'][i][15])
	Stellar_Mass_z5[i] = np.log10(gals_z5['StellarMass'][i]*1.0E10/0.673)

for i in range(0,len(gals_z6['Type'])):
	AGB_Dust_Mass_z6[i] = np.log10(gals_z6['DustMassISM'][i][0]+gals_z6['DustMassISM'][i][1]+gals_z6['DustMassISM'][i][2]+gals_z6['DustMassISM'][i][3])
	SNII_Dust_Mass_z6[i] = np.log10(gals_z6['DustMassISM'][i][4]+gals_z6['DustMassISM'][i][5]+gals_z6['DustMassISM'][i][6]+gals_z6['DustMassISM'][i][7])
	SNIA_Dust_Mass_z6[i] = np.log10(gals_z6['DustMassISM'][i][8]+gals_z6['DustMassISM'][i][9]+gals_z6['DustMassISM'][i][10]+gals_z6['DustMassISM'][i][11])
	Growth_Dust_Mass_z6[i] =np.log10 (gals_z6['DustMassISM'][i][12]+gals_z6['DustMassISM'][i][13]+gals_z6['DustMassISM'][i][14]+gals_z6['DustMassISM'][i][15])
	Stellar_Mass_z6[i] = np.log10(gals_z6['StellarMass'][i]*1.0E10/0.673)

for i in range(0,len(gals_z7['Type'])):
	AGB_Dust_Mass_z7[i] = np.log10(gals_z7['DustMassISM'][i][0]+gals_z7['DustMassISM'][i][1]+gals_z7['DustMassISM'][i][2]+gals_z7['DustMassISM'][i][3])
	SNII_Dust_Mass_z7[i] =np.log10 (gals_z7['DustMassISM'][i][4]+gals_z7['DustMassISM'][i][5]+gals_z7['DustMassISM'][i][6]+gals_z7['DustMassISM'][i][7])
	SNIA_Dust_Mass_z7[i] = np.log10(gals_z7['DustMassISM'][i][8]+gals_z7['DustMassISM'][i][9]+gals_z7['DustMassISM'][i][10]+gals_z7['DustMassISM'][i][11])
	Growth_Dust_Mass_z7[i] =np.log10 (gals_z7['DustMassISM'][i][12]+gals_z7['DustMassISM'][i][13]+gals_z7['DustMassISM'][i][14]+gals_z7['DustMassISM'][i][15])
	Stellar_Mass_z7[i] = np.log10(gals_z7['StellarMass'][i]*1.0E10/0.673)

for i in range(0,len(gals_z8['Type'])):
	AGB_Dust_Mass_z8[i] = np.log10(gals_z8['DustMassISM'][i][0]+gals_z8['DustMassISM'][i][1]+gals_z8['DustMassISM'][i][2]+gals_z8['DustMassISM'][i][3])
	SNII_Dust_Mass_z8[i] =np.log10 (gals_z8['DustMassISM'][i][4]+gals_z8['DustMassISM'][i][5]+gals_z8['DustMassISM'][i][6]+gals_z8['DustMassISM'][i][7])
	SNIA_Dust_Mass_z8[i] = np.log10(gals_z8['DustMassISM'][i][8]+gals_z8['DustMassISM'][i][9]+gals_z8['DustMassISM'][i][10]+gals_z8['DustMassISM'][i][11])
	Growth_Dust_Mass_z8[i] = np.log10(gals_z8['DustMassISM'][i][12]+gals_z8['DustMassISM'][i][13]+gals_z8['DustMassISM'][i][14]+gals_z8['DustMassISM'][i][15])
	Stellar_Mass_z8[i] = np.log10(gals_z8['StellarMass'][i]*1.0E10/0.673)

for i in range(0,len(gals_z9['Type'])):
	AGB_Dust_Mass_z9[i] = np.log10(gals_z9['DustMassISM'][i][0]+gals_z9['DustMassISM'][i][1]+gals_z9['DustMassISM'][i][2]+gals_z9['DustMassISM'][i][3])
	SNII_Dust_Mass_z9[i] = np.log10(gals_z9['DustMassISM'][i][4]+gals_z9['DustMassISM'][i][5]+gals_z9['DustMassISM'][i][6]+gals_z9['DustMassISM'][i][7])
	SNIA_Dust_Mass_z9[i] = np.log10(gals_z9['DustMassISM'][i][8]+gals_z9['DustMassISM'][i][9]+gals_z9['DustMassISM'][i][10]+gals_z9['DustMassISM'][i][11])
	Growth_Dust_Mass_z9[i] = np.log10(gals_z9['DustMassISM'][i][12]+gals_z9['DustMassISM'][i][13]+gals_z9['DustMassISM'][i][14]+gals_z9['DustMassISM'][i][15])
	Stellar_Mass_z9[i] = np.log10(gals_z9['StellarMass'][i]*1.0E10/0.673)



redshift = np.array([0,1,2,3,4,5,6,7,8,9])
avg_AGB=np.zeros(10)
sdv_AGB=np.zeros(10)
avg_SNII=np.zeros(10)
sdv_SNII=np.zeros(10)
avg_GROW=np.zeros(10)
sdv_GROW=np.zeros(10)
avg_SNIA=np.zeros(10)
sdv_SNIA=np.zeros(10)


avg_AGB[0] = (np.median(AGB_Dust_Mass_z0[(Stellar_Mass_z0>=9.0)]))
sdv_AGB[0] = (np.std(AGB_Dust_Mass_z0[(Stellar_Mass_z0>=9.0)]))
avg_AGB[1] = (np.median(AGB_Dust_Mass_z1[(Stellar_Mass_z1>=9.0)]))
sdv_AGB[1] = (np.std(AGB_Dust_Mass_z1[(Stellar_Mass_z1>=9.0)]))
avg_AGB[2] = (np.median(AGB_Dust_Mass_z2[(Stellar_Mass_z2>=9.0)]))
sdv_AGB[2] = (np.std(AGB_Dust_Mass_z2[(Stellar_Mass_z2>=9.0)]))
avg_AGB[3] = (np.median(AGB_Dust_Mass_z3[(Stellar_Mass_z3>=9.0)]))
sdv_AGB[3] = (np.std(AGB_Dust_Mass_z3[(Stellar_Mass_z3>=9.0)]))
avg_AGB[4] = (np.median(AGB_Dust_Mass_z4[(Stellar_Mass_z4>=9.0)]))
sdv_AGB[4] = (np.std(AGB_Dust_Mass_z4[(Stellar_Mass_z4>=9.0)]))
avg_AGB[5] = (np.median(AGB_Dust_Mass_z5[(Stellar_Mass_z5>=9.0)]))
sdv_AGB[5] = (np.std(AGB_Dust_Mass_z5[(Stellar_Mass_z5>=9.0)]))
avg_AGB[6] = (np.median(AGB_Dust_Mass_z6[(Stellar_Mass_z6>=9.0)]))
sdv_AGB[6] = (np.std(AGB_Dust_Mass_z6[(Stellar_Mass_z6>=9.0)]))
avg_AGB[7] = (np.median(AGB_Dust_Mass_z7[(Stellar_Mass_z7>=9.0)]))
sdv_AGB[7] = (np.std(AGB_Dust_Mass_z7[(Stellar_Mass_z7>=9.0)]))
avg_AGB[8] = (np.median(AGB_Dust_Mass_z8[(Stellar_Mass_z8>=9.0)]))
sdv_AGB[8] = (np.std(AGB_Dust_Mass_z8[(Stellar_Mass_z8>=9.0)]))
avg_AGB[9] = (np.median(AGB_Dust_Mass_z9[(Stellar_Mass_z9>=9.0)]))
sdv_AGB[9] = (np.std(AGB_Dust_Mass_z9[(Stellar_Mass_z9>=9.0)]))

avg_SNII[0] = (np.median(SNII_Dust_Mass_z0[(Stellar_Mass_z0>=9.0)]))
sdv_SNII[0] = (np.std(SNII_Dust_Mass_z0[(Stellar_Mass_z0>=9.0)]))
avg_SNII[1] = (np.median(SNII_Dust_Mass_z1[(Stellar_Mass_z1>=9.0)]))
sdv_SNII[1] = (np.std(SNII_Dust_Mass_z1[(Stellar_Mass_z1>=9.0)]))
avg_SNII[2] = (np.median(SNII_Dust_Mass_z2[(Stellar_Mass_z2>=9.0)]))
sdv_SNII[2] = (np.std(SNII_Dust_Mass_z2[(Stellar_Mass_z2>=9.0)]))
avg_SNII[3] = (np.median(SNII_Dust_Mass_z3[(Stellar_Mass_z3>=9.0)]))
sdv_SNII[3] = (np.std(SNII_Dust_Mass_z3[(Stellar_Mass_z3>=9.0)]))
avg_SNII[4] = (np.median(SNII_Dust_Mass_z4[(Stellar_Mass_z4>=9.0)]))
sdv_SNII[4] = (np.std(SNII_Dust_Mass_z4[(Stellar_Mass_z4>=9.0)]))
avg_SNII[5] = (np.median(SNII_Dust_Mass_z5[(Stellar_Mass_z5>=9.0)]))
sdv_SNII[5] = (np.std(SNII_Dust_Mass_z5[(Stellar_Mass_z5>=9.0)]))
avg_SNII[6] = (np.median(SNII_Dust_Mass_z6[(Stellar_Mass_z6>=9.0)]))
sdv_SNII[6] = (np.std(SNII_Dust_Mass_z6[(Stellar_Mass_z6>=9.0)]))
avg_SNII[7] = (np.median(SNII_Dust_Mass_z7[(Stellar_Mass_z7>=9.0)]))
sdv_SNII[7] = (np.std(SNII_Dust_Mass_z7[(Stellar_Mass_z7>=9.0)]))
avg_SNII[8] = (np.median(SNII_Dust_Mass_z8[(Stellar_Mass_z8>=9.0)]))
sdv_SNII[8] = (np.std(SNII_Dust_Mass_z8[(Stellar_Mass_z8>=9.0)]))
avg_SNII[9] = (np.median(SNII_Dust_Mass_z9[(Stellar_Mass_z9>=9.0)]))
sdv_SNII[9] = (np.std(SNII_Dust_Mass_z9[(Stellar_Mass_z9>=9.0)]))

avg_SNIA[0] = (np.median(SNIA_Dust_Mass_z0[(Stellar_Mass_z0>=9.0)]))
sdv_SNIA[0] = (np.std(SNIA_Dust_Mass_z0[(Stellar_Mass_z0>=9.0)]))
avg_SNIA[1] = (np.median(SNIA_Dust_Mass_z1[(Stellar_Mass_z1>=9.0)]))
sdv_SNIA[1] = (np.std(SNIA_Dust_Mass_z1[(Stellar_Mass_z1>=9.0)]))
avg_SNIA[2] = (np.median(SNIA_Dust_Mass_z2[(Stellar_Mass_z2>=9.0)]))
sdv_SNIA[2] = (np.std(SNIA_Dust_Mass_z2[(Stellar_Mass_z2>=9.0)]))
avg_SNIA[3] = (np.median(SNIA_Dust_Mass_z3[(Stellar_Mass_z3>=9.0)]))
sdv_SNIA[3] = (np.std(SNIA_Dust_Mass_z3[(Stellar_Mass_z3>=9.0)]))
avg_SNIA[4] = (np.median(SNIA_Dust_Mass_z4[(Stellar_Mass_z4>=9.0)]))
sdv_SNIA[4] = (np.std(SNIA_Dust_Mass_z4[(Stellar_Mass_z4>=9.0)]))
avg_SNIA[5] = (np.median(SNIA_Dust_Mass_z5[(Stellar_Mass_z5>=9.0)]))
sdv_SNIA[5] = (np.std(SNIA_Dust_Mass_z5[(Stellar_Mass_z5>=9.0)]))
avg_SNIA[6] = (np.median(SNIA_Dust_Mass_z6[(Stellar_Mass_z6>=9.0)]))
sdv_SNIA[6] = (np.std(SNIA_Dust_Mass_z6[(Stellar_Mass_z6>=9.0)]))
avg_SNIA[7] = (np.median(SNIA_Dust_Mass_z7[(Stellar_Mass_z7>=9.0)]))
sdv_SNIA[7] = (np.std(SNIA_Dust_Mass_z7[(Stellar_Mass_z7>=9.0)]))
avg_SNIA[8] = (np.median(SNIA_Dust_Mass_z8[(Stellar_Mass_z8>=9.0)]))
sdv_SNIA[8] = (np.std(SNIA_Dust_Mass_z8[(Stellar_Mass_z8>=9.0)]))
avg_SNIA[9] = (np.median(SNIA_Dust_Mass_z9[(Stellar_Mass_z9>=9.0)]))
sdv_SNIA[9] = (np.std(SNIA_Dust_Mass_z9[(Stellar_Mass_z9>=9.0)]))

avg_GROW[0] = (np.median(Growth_Dust_Mass_z0[(Stellar_Mass_z0>=9.0)]))
sdv_GROW[0] = (np.std(Growth_Dust_Mass_z0[(Stellar_Mass_z0>=9.0)]))
avg_GROW[1] = (np.median(Growth_Dust_Mass_z1[(Stellar_Mass_z1>=9.0)]))
sdv_GROW[1] = (np.std(Growth_Dust_Mass_z1[(Stellar_Mass_z1>=9.0)]))
avg_GROW[2] = (np.median(Growth_Dust_Mass_z2[(Stellar_Mass_z2>=9.0)]))
sdv_GROW[2] = (np.std(Growth_Dust_Mass_z2[(Stellar_Mass_z2>=9.0)]))
avg_GROW[3] = (np.median(Growth_Dust_Mass_z3[(Stellar_Mass_z3>=9.0)]))
sdv_GROW[3] = (np.std(Growth_Dust_Mass_z3[(Stellar_Mass_z3>=9.0)]))
avg_GROW[4] = (np.median(Growth_Dust_Mass_z4[(Stellar_Mass_z4>=9.0)]))
sdv_GROW[4] = (np.std(Growth_Dust_Mass_z4[(Stellar_Mass_z4>=9.0)]))
avg_GROW[5] = (np.median(Growth_Dust_Mass_z5[(Stellar_Mass_z5>=9.0)]))
sdv_GROW[5] = (np.std(Growth_Dust_Mass_z5[(Stellar_Mass_z5>=9.0)]))
avg_GROW[6] = (np.median(Growth_Dust_Mass_z6[(Stellar_Mass_z6>=9.0)]))
sdv_GROW[6] = (np.std(Growth_Dust_Mass_z6[(Stellar_Mass_z6>=9.0)]))
avg_GROW[7] = (np.median(Growth_Dust_Mass_z7[(Stellar_Mass_z7>=9.0)]))
sdv_GROW[7] = (np.std(Growth_Dust_Mass_z7[(Stellar_Mass_z7>=9.0)]))
avg_GROW[8] = (np.median(Growth_Dust_Mass_z8[(Stellar_Mass_z8>=9.0)]))
sdv_GROW[8] = (np.std(Growth_Dust_Mass_z8[(Stellar_Mass_z8>=9.0)]))
avg_GROW[9] = (np.median(Growth_Dust_Mass_z9[(Stellar_Mass_z9>=9.0)]))
sdv_GROW[9] = (np.std(Growth_Dust_Mass_z9[(Stellar_Mass_z9>=9.0)]))




print avg_GROW[5], sdv_GROW[5]


plt.xlabel(r'redshift', fontsize=14,labelpad=10)
plt.ylabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=0)
plt.tick_params(axis='both', which='major', labelsize=10)
plt.tick_params(axis='both', which='minor', labelsize=8)
plt.errorbar(redshift,avg_AGB,sdv_AGB, color='r',label='AGB')
plt.errorbar(redshift,avg_SNII,sdv_SNII, color='b',label='SNII')
plt.errorbar(redshift,avg_GROW,sdv_GROW, color='g',label='Grown')
plt.errorbar(redshift,avg_SNIA,sdv_SNIA, color='k',label='SNIA')
plt.legend(loc='lower left')
pylab.savefig('./graphs/dust_redshift.png', bbox_inches=0)

#plt.show()


