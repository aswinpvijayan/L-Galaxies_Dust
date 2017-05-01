# Script to read in pickled file of galaxy properties.
# Requires pickleFile to hold name of desired input file.
# Returns galaxy properties in gals.

def fit_scatter(x, y, ret_n=False, ret_sterr=False, ret_median=False, nbins=10):
    '''
    Bins scattered points and fits with error bars
    '''
    import numpy as np
    import scipy.stats

    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)

    bin_centres = (_[1:] + _[:-1])/2.

    if ret_sterr:
        stderr = std/np.sqrt(n)
        if ret_n:
            if ret_median:
                median = scipy.stats.binned_statistic(x, y, statistic='median', bins=nbins)[0]
                mederr = std*1.2533
                return bin_centres, mean, std, stderr, n, median, mederr
            return bin_centres, mean, std, stderr, n
        return bin_centres, mean, std, stderr

    if ret_n:
        if ret_median:
            median = scipy.stats.binned_statistic(x, y, statistic='median', bins=nbins)
            mederr = std*1.2533
            return bin_centres, mean, std, n, median, mederr
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
Remy_Dust_Metal_Ratio = Remy_DM / ((10**(Remy_Metals - 8.69 + 0.0134))*(Remy_MHI + Remy_H2mw))


####RR2015 z=0

RR_SM,RR_SMerr,RR_DM1,RR_DM1_up,RR_DM1_down,RR_DM2, RR_DM2_up, RR_DM2_down, RR_MHI, RR_MHI_err, RR_Oxygen, RR_MH2_MW, RR_MH2_Z 	= np.loadtxt('../observations/RR2015.txt',unpack=True,comments='#')

RR_DTG1A = RR_DM1 - np.log10(RR_MHI + RR_MH2_Z)
RR_DTG1B = RR_DM1 - np.log10(RR_MHI + RR_MH2_MW)
RR_DTG2A = RR_DM2 - np.log10(RR_MHI + RR_MH2_Z)
RR_DTG2B = RR_DM2 - np.log10(RR_MHI + RR_MH2_MW)

RR_DTM1A = RR_DM1 - np.log10((10**(RR_Oxygen - 8.69 - 1.87289520164))*(RR_MHI + RR_MH2_Z))
RR_DTM1B = RR_DM1 - np.log10((10**(RR_Oxygen - 8.69 - 1.87289520164))*(RR_MHI + RR_MH2_MW))
RR_DTM2A = RR_DM2 - np.log10((10**(RR_Oxygen - 8.69 - 1.87289520164))*(RR_MHI + RR_MH2_Z))
RR_DTM2B = RR_DM2 - np.log10((10**(RR_Oxygen - 8.69 - 1.87289520164))*(RR_MHI + RR_MH2_MW))

#-1.87289520164
#+ 0.0134
print len(RR_SM),len(RR_DTM1A),len(RR_DTM1B),len(RR_DTM2A),len(RR_DTM2B)

for i in range(0, len(RR_SM)):
	print RR_SM[i],RR_DTM1A[i],RR_DTM1B[i],RR_DTM2A[i],RR_DTM2B[i]

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
#     if Ciesla_ID1[i] != Ciesla_ID2[i]:
#         print Ciesla_ID1[i], Ciesla_ID2[i]

#### Wiseman2016
Wiseman_z_z2, Wiseman_SM_z2, Wiseman_SMuperr_z2, Wiseman_SMdownerr_z2, Wiseman_SFR_z2, Wiseman_SFRerr_z2, Wiseman_Metals_z2, Wiseman_Metalserr_z2, Wiseman_DTM_z2, Wiseman_DTMerr_z2 = np.loadtxt('../observations/wiseman2016_z2.txt',unpack=True,comments='#')
Wiseman_z_z3, Wiseman_SM_z3, Wiseman_SMuperr_z3, Wiseman_SMdownerr_z3, Wiseman_SFR_z3, Wiseman_SFRerr_z3, Wiseman_Metals_z3, Wiseman_Metalserr_z3, Wiseman_DTM_z3, Wiseman_DTMerr_z3 = np.loadtxt('../observations/wiseman2016_z3.txt',unpack=True,comments='#')
Wiseman_z_z4, Wiseman_SM_z4, Wiseman_SMuperr_z4, Wiseman_SMdownerr_z4, Wiseman_SFR_z4, Wiseman_SFRerr_z4, Wiseman_Metals_z4, Wiseman_Metalserr_z4, Wiseman_DTM_z4, Wiseman_DTMerr_z4 = np.loadtxt('../observations/wiseman2016_z4.txt',unpack=True,comments='#')



print "Redshift [Number of galaxies in each mass bin]"    

for loop in range(0,10):
#for loop in range(0,1):

    #Read in L-Galaxies data
    fin = open('../data/'+str(sys.argv[1])+'/lgal_z'+str(loop)+'.pkl','rb')
    gals=cPickle.load(fin)
    fin.close()

    Type = np.zeros(len(gals['Type']))
    #print Type

    Stellar_Mass = np.zeros(len(gals['Type']))
    ColdGas = np.zeros(len(gals['Type']))
    
    SFR = np.zeros(len(gals['Type']))
    sSFR = np.zeros(len(gals['Type']))
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
    
    H_Mass_Dust = np.zeros(len(gals['Type']))
    He_Mass_Dust = np.zeros(len(gals['Type']))
    Cb_Mass_Dust = np.zeros(len(gals['Type']))
    N_Mass_Dust = np.zeros(len(gals['Type']))
    O_Mass_Dust = np.zeros(len(gals['Type']))
    Ne_Mass_Dust = np.zeros(len(gals['Type']))
    Mg_Mass_Dust = np.zeros(len(gals['Type']))
    Si_Mass_Dust = np.zeros(len(gals['Type']))
    S_Mass_Dust = np.zeros(len(gals['Type']))
    Ca_Mass_Dust = np.zeros(len(gals['Type']))
    Fe_Mass_Dust = np.zeros(len(gals['Type']))

    New_Dust_Mass = np.zeros(len(gals['Type']))
    
    DustRate_AGB  = np.zeros(len(gals['Type']))
    DustRate_SNII = np.zeros(len(gals['Type']))
    DustRate_SNIA = np.zeros(len(gals['Type']))
    DustRate_GROW = np.zeros(len(gals['Type']))
    DustRate_DEST = np.zeros(len(gals['Type']))




    for i in range(0,len(gals['Type'])):
    	for j in range(0,11):
            New_Dust_Mass[i] += gals['Dust_elements'][i][j]  
        for j in range(2,11):
            Metals[i] += gals['ColdGas_elements'][i][j]
        
        DustRate_AGB[i]  = gals['DustRatesISM'][i][0]
        DustRate_SNII[i] = gals['DustRatesISM'][i][1]
        DustRate_SNIA[i] = gals['DustRatesISM'][i][2]
        DustRate_GROW[i] = gals['DustRatesISM'][i][3]
        DustRate_DEST[i] = gals['DustRatesISM'][i][4]
        
        Metals_AGB[i]  = gals['MetalsColdGas'][i][0]*1.0E10/0.673
        Metals_SNII[i] = gals['MetalsColdGas'][i][1]*1.0E10/0.673
        Metals_SNIA[i] = gals['MetalsColdGas'][i][2]*1.0E10/0.673
        
        H_Mass[i] = gals['ColdGas_elements'][i][0]
        He_Mass[i]= gals['ColdGas_elements'][i][1]
        Cb_Mass[i]= gals['ColdGas_elements'][i][2]
        N_Mass[i] = gals['ColdGas_elements'][i][3]
        O_Mass[i] = gals['ColdGas_elements'][i][4]
        Ne_Mass[i]= gals['ColdGas_elements'][i][5]
        Mg_Mass[i]= gals['ColdGas_elements'][i][6]
        Si_Mass[i]= gals['ColdGas_elements'][i][7]
        S_Mass[i] = gals['ColdGas_elements'][i][8]
        Ca_Mass[i]= gals['ColdGas_elements'][i][9]
        Fe_Mass[i]= gals['ColdGas_elements'][i][10]
        
        H_Mass_Dust[i] = gals['Dust_elements'][i][0]
        He_Mass_Dust[i]= gals['Dust_elements'][i][1]
        Cb_Mass_Dust[i]= gals['Dust_elements'][i][2]
        N_Mass_Dust[i] = gals['Dust_elements'][i][3]
        O_Mass_Dust[i] = gals['Dust_elements'][i][4]
        Ne_Mass_Dust[i]= gals['Dust_elements'][i][5]
        Mg_Mass_Dust[i]= gals['Dust_elements'][i][6]
        Si_Mass_Dust[i]= gals['Dust_elements'][i][7]
        S_Mass_Dust[i] = gals['Dust_elements'][i][8]
        Ca_Mass_Dust[i]= gals['Dust_elements'][i][9]
        Fe_Mass_Dust[i]= gals['Dust_elements'][i][10]
        
    Stellar_Mass = gals['StellarMass']*1.0E10/0.673
    ColdGas = gals['ColdGas']*1.0E10/0.673
    SFR = gals['Sfr']
    Type = gals['Type']
    sSFR = SFR / Stellar_Mass
        
#    This doesnt work, needs to be done like the dust        
#     H_Mass = gals['ColdGas_elements'][0]
#     He_Mass = gals['ColdGas_elements'][1]
#     Cb_Mass = gals['ColdGas_elements'][2]
#     N_Mass = gals['ColdGas_elements'][3]
#     O_Mass = gals['ColdGas_elements'][4]
#     Ne_Mass = gals['ColdGas_elements'][5]
#     Mg_Mass = gals['ColdGas_elements'][6]
#     Si_Mass = gals['ColdGas_elements'][7]
#     S_Mass = gals['ColdGas_elements'][8]
#     Ca_Mass = gals['ColdGas_elements'][9]
#     Fe_Mass = gals['ColdGas_elements'][10]
            
    #Stellar_Mass_Condition = 1.0E9
    Stellar_Mass_Condition = 0.0
    Dust_Mass_Condition = 0.0 
        
    #--------------------- DATA CHECKS    
    
    if np.ptp(Type, axis=0) != 2:
        print "Min(Type) =",min(Type)
        print "Max(Type) =",max(Type)
        exit("Error with snap template. Type does not equal 0,1,2.")
        
    if min(Stellar_Mass)<0.0 or min(Metals_AGB)<0.0 or min(Metals_SNII)<0.0 or min(Metals_SNIA)<0.0:
        print "Min(SM) =",min(Stellar_Mass)
        print "Min(Metals_AGB) =", min(Metals_AGB)
        print "Min(Metals_SNII) =", min(Metals_SNII)
        print "Min(Metals_SNIA) =", min(Metals_SNIA)
        exit("Negative mass in metallicity or stellarmass detected")
        
    if min(H_Mass) <0.0 or min(He_Mass)<0.0 or min(Cb_Mass)<0.0 or min(N_Mass)<0.0 or min(O_Mass)<0.0 or min(Ne_Mass)<0.0 or min(Mg_Mass)<0.0 or min(Si_Mass)<0.0 or min(S_Mass)<0.0 or min(Ca_Mass)<0.0 or min(Fe_Mass)<0.0:
        print "Min(H) =",min(H_Mass)
        print "Min(He) =",min(He_Mass)
        print "Min(Cb) =",min(Cb_Mass)
        print "Min(N) =",min(N_Mass)
        print "Min(O) =",min(O_Mass)
        print "Min(Ne) =",min(Ne_Mass)
        print "Min(Mg) =",min(Mg_Mass)
        print "Min(Si) =",min(Si_Mass)
        print "Min(S) =",min(S_Mass)
        print "Min(Ca) =",min(Ca_Mass)
        print "Min(Fe) =",min(Fe_Mass)
        exit("Negative mass in metals detected")
    
    if min(H_Mass_Dust) <0.0 or min(He_Mass_Dust)<0.0 or min(Cb_Mass_Dust)<0.0 or min(N_Mass_Dust)<0.0 or min(O_Mass_Dust)<0.0 or min(Ne_Mass_Dust)<0.0 or min(Mg_Mass_Dust)<0.0 or min(Si_Mass_Dust)<0.0 or min(S_Mass_Dust)<0.0 or min(Ca_Mass_Dust)<0.0 or min(Fe_Mass_Dust)<0.0:
        print "Min(H_Mass_Dust) =",min(H_Mass_Dust)
        print "Min(He_Mass_Dust) =",min(He_Mass_Dust)
        print "Min(Cb_Mass_Dust) =",min(Cb_Mass_Dust)
        print "Min(N_Mass_Dust) =",min(N_Mass_Dust)
        print "Min(O_Mass_Dust) =",min(O_Mass_Dust)
        print "Min(Ne_Mass_Dust) =",min(Ne_Mass_Dust)
        print "Min(Mg_Mass_Dust) =",min(Mg_Mass_Dust)
        print "Min(Si_Mass_Dust) =",min(Si_Mass_Dust)
        print "Min(S_Mass_Dust) =",min(S_Mass_Dust)
        print "Min(Ca_Mass_Dust) =",min(Ca_Mass_Dust)
        print "Min(Fe_Mass_Dust) =",min(Fe_Mass_Dust)
        exit("Negative mass in dust detected")
            
    #---------------------All NEW dust Plots

    #condition = np.logical_and(New_Dust_Mass>0,Stellar_Mass>Stellar_Mass_Condition)
    #condition = np.logical_and(np.logical_and(New_Dust_Mass>0, np.log10(SFR) > 0.0),Stellar_Mass>Stellar_Mass_Condition)
    condition = np.logical_and(sSFR > 0.0345E-9, np.logical_and(New_Dust_Mass > 0, Type == 0))
    
    log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
    log_New_Dust_Mass = np.log10(New_Dust_Mass[condition==1])
    
    SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count,Dust_median,Dust_mederr = fit_scatter(log_Stellar_Mass, log_New_Dust_Mass, ret_n=True, ret_sterr=True,ret_median=True, nbins=20)
    print loop, count
    print SM_bins
    print Dust_bins
    print Dust_median
    np.savetxt('./binned_data/SM_DM_z'+str(loop)+'_'+str(sys.argv[1])+'.txt',np.c_[SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count])

    if(sum(count)>0):    
        # You typically want your plot to be ~1.33x wider than tall. This plot is a rare    
        # exception because of the number of lines being plotted on it.    
        # Common sizes: (10, 7.5) and (12, 9)    
        #plt.figure(figsize=(12, 9))            
        fig = plt.figure(figsize=(9,6))
        plt.xlim([6,12])
        plt.ylim([0,10.2])
        plt.hexbin(log_Stellar_Mass,log_New_Dust_Mass,gridsize=500,mincnt=1, label='L-Galaxies 2Dhist')
        plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
        plt.ylabel(r'log$_{10}$(M$_{\rm{d}}$/M$_{\odot}$)', fontsize=18)
        plt.tick_params(axis='both', which='major', labelsize=12,width=2)
        plt.tick_params(axis='both', which='minor', labelsize=12,width=2)
        #plt.text(7,1,"N = "+str(sum(count))+"\nz = "+str(loop))
        plt.text(9,1,"z = "+str(loop)+"\nN = "+str(sum(count)),fontsize=16)
        if(loop == 0):
            plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='g',label='RemyRuyer2014',fmt='o')
            #plt.errorbar(np.log10(Bourne_MEDMSTAR), np.log10(Bourne_MEDMDUST), yerr = np.log10(Bourne_MEDMDUSTERR/Bourne_MEDMDUST) , color='orange',label='Bourne2012',fmt='o')
            #plt.errorbar(Ciesla_SM, Ciesla_DM, yerr = Ciesla_DMerr , color='r',label='Ciesla2014',fmt='o')
            #plt.errorbar(Santini_SM_z0, Santini_DM_z0, yerr = (Santini_DMdownerr_z0, Santini_DMuperr_z0), color='b',label='Santini2014',fmt='o')
            plt.errorbar(RR_SM, RR_DM1, yerr = (RR_DM1_down, RR_DM1_up),color='b',label='RR2015-1',fmt='o')
            plt.errorbar(RR_SM, RR_DM2, yerr = (RR_DM2_down, RR_DM2_up),color='r',label='RR2015-2',fmt='o')
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
    
        plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='k',label='L-Galaxies Mean',linewidth=2)
        plt.errorbar(SM_bins,Dust_median,yerr=(Dust_mederr),color='b',label='L-Galaxies Median',linewidth=2)
        
        plt.legend(loc='lower right')
        
        axes = fig.get_axes()
        for ax in axes:
            [i.set_linewidth(2.1) for i in ax.spines.itervalues()]
        
        
        pylab.savefig('./graphs/stellar_Newdust_z'+str(loop)+'_'+str(sys.argv[1])+'.png', bbox_inches=0)
        plt.close()
    
    if loop == 0:
        avg_NEW_dust = np.array([np.mean(log_New_Dust_Mass)])
        std_NEW_dust = np.array([np.std(log_New_Dust_Mass)/np.sqrt(len(log_New_Dust_Mass))])

    else:
        avg_NEW_dust = np.append(avg_NEW_dust,np.mean(log_New_Dust_Mass))
        std_NEW_dust = np.append(std_NEW_dust,np.std(log_New_Dust_Mass)/np.sqrt(len(log_New_Dust_Mass)))
        
        
        
        
        
        

    #---------------------Metal Dust ratios

    #condition = np.logical_and(np.logical_and(Metals>0,New_Dust_Mass>0),Stellar_Mass>0)
    
    log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
    log_New_Dust_Mass = np.log10(New_Dust_Mass[condition==1])
    log_Metals = np.log10(Metals[condition==1])
    
    Ratio = log_New_Dust_Mass - log_Metals
    
    SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count = fit_scatter(log_Stellar_Mass, Ratio, ret_n=True, ret_sterr=True, nbins=10)
    print loop, count

    if(sum(count)>0):    
        plt.xlim([6,12])
        plt.ylim([-3,2])
        plt.hexbin(log_Stellar_Mass,Ratio,gridsize=500,mincnt=1, label='Dust All')
        plt.errorbar(SM_bins,Dust_bins,yerr=(Dust_std_err),color='r',label='Dust/Metal ratio')
        plt.xlabel(r'log$_{10}$(M_d/M$_{\odot}$)', fontsize=14)
        plt.ylabel(r'log$_{10}$(M_d/M_m)', fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.tick_params(axis='both', which='minor', labelsize=8)
        #plt.legend(loc='lower right')
        plt.text(2,-7,"N = "+str(sum(count)))
        plt.text(2,-5,"z = "+str(loop)+" :All dust")
        
        if(loop == 0):
            plt.scatter(Remy_SM, np.log10(Remy_Dust_Metal_Ratio),color='g',label='RR2014',marker='o')
            plt.scatter(RR_SM, (RR_DTM1A), color='b',label='1A',marker='o')
            plt.scatter(RR_SM, (RR_DTM1B), color='cyan',label='1B',marker='o')
            plt.scatter(RR_SM, (RR_DTM2A), color='r',label='2A',marker='o')
            plt.scatter(RR_SM, (RR_DTM2B), color='orange',label='2B',marker='o')
            
            plt.legend(loc='upper right')
        
        pylab.savefig('./graphs/stellar_dustmetalratio_z'+str(loop)+'_'+str(sys.argv[1])+'.png', bbox_inches=0)
        plt.close()


    #---------------------Dust Gas ratios

    #condition = np.logical_and(np.logical_and(ColdGas>0,New_Dust_Mass>0),Stellar_Mass>0)
    
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
        plt.xlabel(r'log$_{10}$(M_d/M$_{\odot}$)', fontsize=14)
        plt.ylabel(r'log$_{10}$(M_d/M_{cg})', fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.tick_params(axis='both', which='minor', labelsize=8)
        #plt.legend(loc='lower right')
        plt.text(2,-7,"N = "+str(sum(count)))
        plt.text(2,-5,"z = "+str(loop)+" :All dust")
        if(loop == 0):
            #plt.errorbar(Remy_SM,np.log10(Remy_DM),yerr=np.log10(Remy_DM)*(Remy_DMerr/100.0),color='g',label='RemyRuyer2014',fmt='o')
            plt.scatter(Remy_SM, np.log10(Remy_Dust_Gas_Ratio),color='g',label='RR2014',marker='o')
            plt.scatter(RR_SM, RR_DTG1A, color='b',label='1A',marker='o') 
            plt.scatter(RR_SM, RR_DTG1B, color='cyan',label='1B',marker='o') 
            plt.scatter(RR_SM, RR_DTG2A, color='r',label='2A',marker='o') 
            plt.scatter(RR_SM, RR_DTG2B, color='orange',label='2B',marker='o') 

            plt.legend(loc='lower right')
        pylab.savefig('./graphs/stellar_dustgasratio_z'+str(loop)+'_'+str(sys.argv[1])+'.png', bbox_inches=0)
        plt.close()
        
        
        
    #---------------------Dust rates

    condition = np.logical_and(np.logical_and(sSFR > 0.0345E-9,DustRate_AGB>0.0),Type == 0)
    log_Stellar_Mass_AGB  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_AGB  = np.log10(DustRate_AGB[condition==1])
    SM_bins_AGB,Dust_bins_AGB,Dust_std_dev_AGB,Dust_std_err_AGB,count_AGB = fit_scatter(log_Stellar_Mass_AGB, log_DustRate_AGB, ret_n=True, ret_sterr=True, nbins=10)
    np.savetxt('./binned_data/dustrates_AGB_z'+str(loop)+'_'+str(sys.argv[1])+'.txt',np.c_[SM_bins_AGB,Dust_bins_AGB,Dust_std_dev_AGB,Dust_std_err_AGB,count_AGB])

    condition = np.logical_and(np.logical_and(sSFR > 0.0345E-9,DustRate_SNII>0.0),Type == 0)
    log_Stellar_Mass_SNII  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_SNII  = np.log10(DustRate_SNII[condition==1])
    SM_bins_SNII,Dust_bins_SNII,Dust_std_dev_SNII,Dust_std_err_SNII,count_SNII = fit_scatter(log_Stellar_Mass_SNII, log_DustRate_SNII, ret_n=True, ret_sterr=True, nbins=10)
    np.savetxt('./binned_data/dustrates_SNII_z'+str(loop)+'_'+str(sys.argv[1])+'.txt',np.c_[SM_bins_SNII,Dust_bins_SNII,Dust_std_dev_SNII,Dust_std_err_SNII,count_SNII])

    condition = np.logical_and(np.logical_and(sSFR > 0.0345E-9,DustRate_SNIA>0.0),Type == 0)
    log_Stellar_Mass_SNIA  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_SNIA  = np.log10(DustRate_SNIA[condition==1])
    SM_bins_SNIA,Dust_bins_SNIA,Dust_std_dev_SNIA,Dust_std_err_SNIA,count_SNIA = fit_scatter(log_Stellar_Mass_SNIA, log_DustRate_SNIA, ret_n=True, ret_sterr=True, nbins=10)
    np.savetxt('./binned_data/dustrates_SNIA_z'+str(loop)+'_'+str(sys.argv[1])+'.txt',np.c_[SM_bins_SNIA,Dust_bins_SNIA,Dust_std_dev_SNIA,Dust_std_err_SNIA,count_SNIA])

    condition = np.logical_and(np.logical_and(sSFR > 0.0345E-9,DustRate_GROW>0.0),Type == 0)
    log_Stellar_Mass_GROW  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_GROW  = np.log10(DustRate_GROW[condition==1])
    SM_bins_GROW,Dust_bins_GROW,Dust_std_dev_GROW,Dust_std_err_GROW,count_GROW = fit_scatter(log_Stellar_Mass_GROW, log_DustRate_GROW, ret_n=True, ret_sterr=True, nbins=10)
    np.savetxt('./binned_data/dustrates_GROW_z'+str(loop)+'_'+str(sys.argv[1])+'.txt',np.c_[SM_bins,Dust_bins,Dust_std_dev,Dust_std_err,count])

    condition = np.logical_and(np.logical_and(sSFR > 0.0345E-9,DustRate_DEST>0.0),Type == 0)
    log_Stellar_Mass_DEST  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_DEST  = np.log10(DustRate_DEST[condition==1])
    SM_bins_DEST,Dust_bins_DEST,Dust_std_dev_DEST,Dust_std_err_DEST,count_DEST = fit_scatter(log_Stellar_Mass_DEST, log_DustRate_DEST, ret_n=True, ret_sterr=True, nbins=10)
    np.savetxt('./binned_data/dustrates_DEST_z'+str(loop)+'_'+str(sys.argv[1])+'.txt',np.c_[SM_bins_DEST,Dust_bins_DEST,Dust_std_dev_DEST,Dust_std_err_DEST,count_DEST])

    
    
    #if(sum(count)>0):    
    plt.xlim([6,12])
    plt.ylim([-10,3])
    #plt.hexbin(log_Stellar_Mass,log_DustRate_AGB,gridsize=500,mincnt=1, label='Dust All')
    plt.errorbar(SM_bins_AGB,Dust_bins_AGB,yerr=(Dust_std_err_AGB),color='b',label='AGB')
    plt.errorbar(SM_bins_SNII,Dust_bins_SNII,yerr=(Dust_std_err_SNII),color='r',label='SNII')
    plt.errorbar(SM_bins_SNIA,Dust_bins_SNIA,yerr=(Dust_std_err_SNIA),color='y',label='SNIA')
    plt.errorbar(SM_bins_GROW,Dust_bins_GROW,yerr=(Dust_std_err_GROW),color='g',label='GROW')
    plt.errorbar(SM_bins_DEST,Dust_bins_DEST,yerr=(Dust_std_err_DEST),color='k',label='DEST')
    
    
    
    plt.xlabel(r'log$_{10}$(M_*/M$_{\odot}$)', fontsize=14,labelpad=10)
    plt.ylabel(r'Rate (M$_{\odot}$/$yr^{-1})$', fontsize=14,labelpad=0)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.tick_params(axis='both', which='minor', labelsize=8)
    plt.legend(loc='lower right')
    #plt.text(2,-7,"N = "+str(sum(count)))
    #plt.text(2,-5,"z = "+str(loop)+" :All dust")
    pylab.savefig('./graphs/dustrates_z'+str(loop)+'_'+str(sys.argv[1])+'.png', bbox_inches=0)
    plt.close()



    #---------------------12+log(O/H) vs. Dust/Metal ratio
    
    
    condition = np.logical_and(np.logical_and(sSFR > 0.0345E-9,O_Mass>0.0), np.logical_and(New_Dust_Mass > 0, Type == 0))

    log_Stellar_Mass = np.log10(Stellar_Mass[condition==1])
    log_New_Dust_Mass = np.log10(New_Dust_Mass[condition==1])
    log_Metals = np.log10(Metals[condition==1])
    Ratio = log_New_Dust_Mass - log_Metals
    
    O_metals = O_Mass[condition==1]
    H_metals = H_Mass[condition==1]
    
    log_metallicity = np.log10( (O_metals/H_metals) * (1.0/16.0) )+ 12.0

    
    Metallicity_bins,Ratio_bins,Ratio_std_dev,Ratio_std_err,count = fit_scatter(log_metallicity, Ratio, ret_n=True, ret_sterr=True, nbins=10)


    if(sum(count)>0):    
        plt.xlim([6,12])
        plt.ylim([-2,1])
        plt.hexbin(log_metallicity,Ratio,gridsize=500,mincnt=1)
        plt.errorbar(Metallicity_bins,Ratio_bins,yerr=(Ratio_std_err),color='r',label='L-Galaxies Mean')        
        plt.xlabel(r'12 + log$_{10}$(O/H)', fontsize=14)
        plt.ylabel(r'log$_{10}$(M$_{d}$/M$_{m}$)', fontsize=14)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.tick_params(axis='both', which='minor', labelsize=8)
        #plt.text(6.2,-1,"N = "+str(sum(count)))
        #plt.text(10,2,"z = "+str(loop)+" :New dust")
        if(loop == 0):
            plt.scatter(Remy_Metals,Remy_Dust_Metal_Ratio,color='g',label='RemyRuyer2014',marker='o')
        if(loop == 2):
            plt.errorbar(Wiseman_Metals_z2, Wiseman_DTM_z2*0.32, yerr = Wiseman_DTMerr_z2, color='g',label='Wiseman2016',fmt='o')
        if(loop == 3):
            plt.errorbar(Wiseman_Metals_z3, Wiseman_DTM_z3*0.32, yerr = Wiseman_DTMerr_z3, color='g',label='Wiseman2016',fmt='o')
        if(loop == 4):
            plt.errorbar(Wiseman_Metals_z4, Wiseman_DTM_z4*0.32, yerr = Wiseman_DTMerr_z4, color='g',label='Wiseman2016',fmt='o')
        plt.legend(loc='lower right')
        
        

        pylab.savefig('./graphs/dustmetalratio_metallicity_z'+str(loop)+'_'+str(sys.argv[1])+'.png', bbox_inches=0)
        plt.close()
        






















