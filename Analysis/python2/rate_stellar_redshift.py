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



#### Read in observational data

#### Remy-Ruyer 2014 z=0
#z0_obs_DM, obs_DM_err, obs_SM = np.loadtxt('../observations/z0.txt',unpack=True)
#z0_obs_DM_act_err = obs_DM*(obs_DM_err/100.0)
#Name DustMass DMerror% SMass MHI MHIerror% 12+log(O/H) MH2,mw MH2,Z
Remy_DM, Remy_DMerr, Remy_SM, Remy_MHI, Remy_MHIerr, Remy_Metals, Remy_H2mw, Remy_MH2z = np.loadtxt('../observations/Remy_Ruyer_2014_KINGFISH_z0.txt',unpack=True,comments='#')
Remy_DM_err_actual = Remy_DM*(Remy_DMerr/100.0)
Remy_Dust_Gas_Ratio = Remy_DM / (Remy_MHI + Remy_H2mw)
Remy_Dust_Metal_Ratio = Remy_DM / (((10**(Remy_Metals - 8.69))*0.0134)*(Remy_MHI + Remy_H2mw))

##print Remy_SM, np.log10(Remy_Dust_Metal_Ratio), Remy_DM, (((10**(Remy_Metals - 8.69))*0.02)*(Remy_MHI + Remy_H2mw))


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
#         #print Ciesla_ID1[i], Ciesla_ID2[i]

#### Wiseman2016
Wiseman_z_z2, Wiseman_SM_z2, Wiseman_SMuperr_z2, Wiseman_SMdownerr_z2, Wiseman_SFR_z2, Wiseman_SFRerr_z2, Wiseman_Metals_z2, Wiseman_Metalserr_z2, Wiseman_DTM_z2, Wiseman_DTMerr_z2 = np.loadtxt('../observations/wiseman2016_z2.txt',unpack=True,comments='#')
Wiseman_z_z3, Wiseman_SM_z3, Wiseman_SMuperr_z3, Wiseman_SMdownerr_z3, Wiseman_SFR_z3, Wiseman_SFRerr_z3, Wiseman_Metals_z3, Wiseman_Metalserr_z3, Wiseman_DTM_z3, Wiseman_DTMerr_z3 = np.loadtxt('../observations/wiseman2016_z3.txt',unpack=True,comments='#')
Wiseman_z_z4, Wiseman_SM_z4, Wiseman_SMuperr_z4, Wiseman_SMdownerr_z4, Wiseman_SFR_z4, Wiseman_SFRerr_z4, Wiseman_Metals_z4, Wiseman_Metalserr_z4, Wiseman_DTM_z4, Wiseman_DTMerr_z4 = np.loadtxt('../observations/wiseman2016_z4.txt',unpack=True,comments='#')

DustRate_z      = np.zeros(9)
DustRate_AGB_z  = np.zeros(9)
DustRate_SNII_z = np.zeros(9)
DustRate_SNIA_z = np.zeros(9)
DustRate_GROW_z = np.zeros(9)
DustRate_DEST_z = np.zeros(9)
DustRate_ALL_z  = np.zeros(9)
SFR_z  = np.zeros(9)

#print "Redshift [Number of galaxies in each mass bin]"    

for loop in range(0,9):
#for loop in range(0,1):
    fin = open('../data/'+str(sys.argv[1])+'/lgal_z'+str(loop)+'.pkl','rb')
    gals=cPickle.load(fin)
    fin.close()
      
    ##############

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

    Dust_Mass = np.zeros(len(gals['Type']))
    
    DustRate_AGB  = np.zeros(len(gals['Type']))
    DustRate_SNII = np.zeros(len(gals['Type']))
    DustRate_SNIA = np.zeros(len(gals['Type']))
    DustRate_GROW = np.zeros(len(gals['Type']))
    DustRate_DEST = np.zeros(len(gals['Type']))
    DustRate_ALL = np.zeros(len(gals['Type']))

    #print "Allocated memory for arrays"
    
    
    Stellar_Mass = gals['StellarMass']*1.0E10/0.673
    ColdGas = gals['ColdGas']*1.0E10/0.673
    SFR = gals['Sfr']
    Type = gals['Type']
    sSFR = SFR / Stellar_Mass

    Dust_Mass = np.sum(gals['Dust_elements'], axis=1)
    
    Metals= (gals['ColdGas_elements'][:,2:11])
    Metal_Mass = np.sum(Metals,axis=1)
    
    #print gals['ColdGas_elements']
    #print Metal_Mass
    #print len(Metal_Mass)
    
    DustRate_AGB  = gals['DustRatesISM'][:,0]
    DustRate_SNII = gals['DustRatesISM'][:,1]
    DustRate_SNIA = gals['DustRatesISM'][:,2]
    DustRate_GROW = gals['DustRatesISM'][:,3]
    DustRate_DEST = gals['DustRatesISM'][:,4]
    DustRate_ALL1 = gals['DustRatesISM'][:,0:4]
    DustRate_ALL  = np.sum(DustRate_ALL1, axis=1)
    
    
    #print gals['DustRatesISM'][0]
    #print DustRate_AGB[0]
    #print DustRate_SNII[0]
    #print DustRate_SNIA[0]
    #print DustRate_GROW[0]
    #print DustRate_DEST[0]

        
    H_Mass  = gals['ColdGas_elements'][:,0]
    He_Mass = gals['ColdGas_elements'][:,1]
    Cb_Mass = gals['ColdGas_elements'][:,2]
    N_Mass  = gals['ColdGas_elements'][:,3]
    O_Mass  = gals['ColdGas_elements'][:,4]
    Ne_Mass = gals['ColdGas_elements'][:,5]
    Mg_Mass = gals['ColdGas_elements'][:,6]
    Si_Mass = gals['ColdGas_elements'][:,7]
    S_Mass  = gals['ColdGas_elements'][:,8]
    Ca_Mass = gals['ColdGas_elements'][:,9]
    Fe_Mass = gals['ColdGas_elements'][:,10]
    
    H_Mass_Dust  = gals['Dust_elements'][:,0]
    He_Mass_Dust = gals['Dust_elements'][:,1]
    Cb_Mass_Dust = gals['Dust_elements'][:,2]
    N_Mass_Dust  = gals['Dust_elements'][:,3]
    O_Mass_Dust  = gals['Dust_elements'][:,4]
    Ne_Mass_Dust = gals['Dust_elements'][:,5]
    Mg_Mass_Dust = gals['Dust_elements'][:,6]
    Si_Mass_Dust = gals['Dust_elements'][:,7]
    S_Mass_Dust  = gals['Dust_elements'][:,8]
    Ca_Mass_Dust = gals['Dust_elements'][:,9]
    Fe_Mass_Dust = gals['Dust_elements'][:,10]
        
    #Stellar_Mass_Condition = 1.0E9
    Stellar_Mass_Condition = 0.0
    Dust_Mass_Condition = 0.0 
        
    #print "Filled arrays with properies"
    #--------------------- DATA CHECKS    
    #print "Data checks"
    
    if np.ptp(Type, axis=0) != 2:
        #print "Min(Type) =",min(Type)
        #print "Max(Type) =",max(Type)
        exit("Error with snap template. Type does not equal 0,1,2.")
        
    if min(Stellar_Mass)<0.0 or min(Metals_AGB)<0.0 or min(Metals_SNII)<0.0 or min(Metals_SNIA)<0.0:
        #print "Min(SM) =",min(Stellar_Mass)
        #print "Min(Metals_AGB) =", min(Metals_AGB)
        #print "Min(Metals_SNII) =", min(Metals_SNII)
        #print "Min(Metals_SNIA) =", min(Metals_SNIA)
        exit("Negative mass in metallicity or stellarmass detected")
        
    if min(H_Mass) <0.0 or min(He_Mass)<0.0 or min(Cb_Mass)<0.0 or min(N_Mass)<0.0 or min(O_Mass)<0.0 or min(Ne_Mass)<0.0 or min(Mg_Mass)<0.0 or min(Si_Mass)<0.0 or min(S_Mass)<0.0 or min(Ca_Mass)<0.0 or min(Fe_Mass)<0.0:
        #print "Min(H) =",min(H_Mass)
        #print "Min(He) =",min(He_Mass)
        #print "Min(Cb) =",min(Cb_Mass)
        #print "Min(N) =",min(N_Mass)
        #print "Min(O) =",min(O_Mass)
        #print "Min(Ne) =",min(Ne_Mass)
        #print "Min(Mg) =",min(Mg_Mass)
        #print "Min(Si) =",min(Si_Mass)
        #print "Min(S) =",min(S_Mass)
        #print "Min(Ca) =",min(Ca_Mass)
        #print "Min(Fe) =",min(Fe_Mass)
        exit("Negative mass in metals detected")
    
    if min(H_Mass_Dust) <0.0 or min(He_Mass_Dust)<0.0 or min(Cb_Mass_Dust)<0.0 or min(N_Mass_Dust)<0.0 or min(O_Mass_Dust)<0.0 or min(Ne_Mass_Dust)<0.0 or min(Mg_Mass_Dust)<0.0 or min(Si_Mass_Dust)<0.0 or min(S_Mass_Dust)<0.0 or min(Ca_Mass_Dust)<0.0 or min(Fe_Mass_Dust)<0.0:
        #print "Min(H_Mass_Dust) =",min(H_Mass_Dust)
        #print "Min(He_Mass_Dust) =",min(He_Mass_Dust)
        #print "Min(Cb_Mass_Dust) =",min(Cb_Mass_Dust)
        #print "Min(N_Mass_Dust) =",min(N_Mass_Dust)
        #print "Min(O_Mass_Dust) =",min(O_Mass_Dust)
        #print "Min(Ne_Mass_Dust) =",min(Ne_Mass_Dust)
        #print "Min(Mg_Mass_Dust) =",min(Mg_Mass_Dust)
        #print "Min(Si_Mass_Dust) =",min(Si_Mass_Dust)
        #print "Min(S_Mass_Dust) =",min(S_Mass_Dust)
        #print "Min(Ca_Mass_Dust) =",min(Ca_Mass_Dust)
        #print "Min(Fe_Mass_Dust) =",min(Fe_Mass_Dust)
        exit("Negative mass in dust detected")
            
    #print "Data checks complete"


    #---------------------dustmetalratio_oxygen Plots
    
    #0.0345E-9
    condition = np.logical_and(np.logical_and(sSFR > 0.0,DustRate_AGB>0.0),Type == 0)
    log_Stellar_Mass_AGB  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_AGB  = np.log10(DustRate_AGB[condition==1])
    SM_bins_AGB,Dust_bins_AGB,Dust_std_dev_AGB,Dust_std_err_AGB,count_AGB, Dust_med_AGB, Dust_mederr_AGB = fit_scatter(log_Stellar_Mass_AGB, log_DustRate_AGB, ret_n=True, ret_sterr=True, ret_median=True, nbins=20)
    #np.savetxt('./binned/dustrates_AGB_z'+str(loop)+'_ALL.txt',np.c_[SM_bins_AGB,Dust_bins_AGB,Dust_std_dev_AGB,Dust_std_err_AGB,count_AGB])

    condition = np.logical_and(np.logical_and(sSFR > 0.0,DustRate_SNII>0.0),Type == 0)
    log_Stellar_Mass_SNII  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_SNII  = np.log10(DustRate_SNII[condition==1])
    SM_bins_SNII,Dust_bins_SNII,Dust_std_dev_SNII,Dust_std_err_SNII,count_SNII, Dust_med_SNII, Dust_mederr_SNII = fit_scatter(log_Stellar_Mass_SNII, log_DustRate_SNII, ret_n=True, ret_sterr=True, ret_median=True,nbins=20)
    #np.savetxt('./binned/dustrates_SNII_z'+str(loop)+'_ALL.txt',np.c_[SM_bins_SNII,Dust_bins_SNII,Dust_std_dev_SNII,Dust_std_err_SNII,count_SNII])

    condition = np.logical_and(np.logical_and(sSFR > 0.0,DustRate_SNIA>0.0),Type == 0)
    log_Stellar_Mass_SNIA  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_SNIA  = np.log10(DustRate_SNIA[condition==1])
    SM_bins_SNIA,Dust_bins_SNIA,Dust_std_dev_SNIA,Dust_std_err_SNIA,count_SNIA, Dust_med_SNIA, Dust_mederr_SNIA = fit_scatter(log_Stellar_Mass_SNIA, log_DustRate_SNIA, ret_n=True, ret_sterr=True,ret_median=True, nbins=20)
    #np.savetxt('./binned/dustrates_SNIA_z'+str(loop)+'_ALL.txt',np.c_[SM_bins_SNIA,Dust_bins_SNIA,Dust_std_dev_SNIA,Dust_std_err_SNIA,count_SNIA])

    condition = np.logical_and(np.logical_and(sSFR > 0.0,DustRate_GROW>0.0),Type == 0)
    log_Stellar_Mass_GROW  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_GROW  = np.log10(DustRate_GROW[condition==1])
    SM_bins_GROW,Dust_bins_GROW,Dust_std_dev_GROW,Dust_std_err_GROW,count_GROW, Dust_med_GROW, Dust_mederr_GROW = fit_scatter(log_Stellar_Mass_GROW, log_DustRate_GROW, ret_n=True, ret_sterr=True,ret_median=True, nbins=20)
    #np.savetxt('./binned/dustrates_GROW_z'+str(loop)+'_ALL.txt',np.c_[SM_bins_GROW,Dust_bins_GROW,Dust_std_dev_GROW,Dust_std_err_GROW,count_GROW])

    condition = np.logical_and(np.logical_and(sSFR > 0.0,DustRate_DEST>0.0),Type == 0)
    log_Stellar_Mass_DEST  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_DEST  = np.log10(DustRate_DEST[condition==1])
    SM_bins_DEST,Dust_bins_DEST,Dust_std_dev_DEST,Dust_std_err_DEST,count_DEST, Dust_med_DEST, Dust_mederr_DEST = fit_scatter(log_Stellar_Mass_DEST, log_DustRate_DEST, ret_n=True, ret_sterr=True,ret_median=True, nbins=20)
    #np.savetxt('./binned/dustrates_DEST_z'+str(loop)+'_ALL.txt',np.c_[SM_bins_DEST,Dust_bins_DEST,Dust_std_dev_DEST,Dust_std_err_DEST,count_DEST])

    condition = np.logical_and(np.logical_and(sSFR > 0.0,DustRate_ALL>0.0),Type == 0)
    log_Stellar_Mass_ALL  = np.log10(Stellar_Mass[condition==1])
    log_DustRate_ALL  = np.log10(DustRate_ALL[condition==1])
    SM_bins_ALL,Dust_bins_ALL,Dust_std_dev_ALL,Dust_std_err_ALL,count_ALL, Dust_med_ALL, Dust_mederr_ALL = fit_scatter(log_Stellar_Mass_ALL, log_DustRate_ALL, ret_n=True, ret_sterr=True,ret_median=True, nbins=20)
    #np.savetxt('./binned/dustrates_ALL_z'+str(loop)+'_ALL.txt',np.c_[SM_bins_ALL,Dust_bins_ALL,Dust_std_dev_ALL,Dust_std_err_ALL,count_ALL])



#     DustRate_z[loop] = loop
#     condition = np.logical_and(Stellar_Mass > (9.0),DustRate_AGB>0.0)
#     DustRate_AGB_z[loop] = np.median(np.log10(DustRate_AGB[condition==1]))
#     condition = np.logical_and(Stellar_Mass > (9.0),DustRate_SNII>0.0)
#     DustRate_SNII_z[loop]= np.median(np.log10(DustRate_SNII[condition==1]))
#     condition = np.logical_and(Stellar_Mass > (9.0),DustRate_SNIA>0.0)
#     DustRate_SNIA_z[loop]= np.median(np.log10(DustRate_SNIA[condition==1]))
#     condition = np.logical_and(Stellar_Mass > (9.0),DustRate_GROW>0.0)
#     DustRate_GROW_z[loop]= np.median(np.log10(DustRate_GROW[condition==1]))
#     condition = np.logical_and(Stellar_Mass > (9.0),DustRate_DEST>0.0)
#     DustRate_DEST_z[loop]= np.median(np.log10(DustRate_DEST[condition==1]))
#     condition = np.logical_and(Stellar_Mass > (9.0),DustRate_ALL>0.0)
#     DustRate_ALL_z[loop] = np.median(np.log10(DustRate_ALL[condition==1]))
#     condition = np.logical_and(Stellar_Mass > (9.0),DustRate_ALL>0.0)
#     SFR_z[loop] = np.median(np.log10(SFR[condition==1]))
    volume = ((480.279/512)/0.673)**3.0
    #volume = (96.1044/0.673)**3.0
    DustRate_z[loop] = loop
    condition = np.logical_and(Stellar_Mass > (9.0),DustRate_AGB>0.0)
    DustRate_AGB_z[loop] = np.log10(np.sum(DustRate_AGB[condition==1])/volume)
    condition = np.logical_and(Stellar_Mass > (9.0),DustRate_SNII>0.0)
    DustRate_SNII_z[loop]= np.log10(np.sum(DustRate_SNII[condition==1])/volume)
    condition = np.logical_and(Stellar_Mass > (9.0),DustRate_SNIA>0.0)
    DustRate_SNIA_z[loop]= np.log10(np.sum(DustRate_SNIA[condition==1])/volume)
    condition = np.logical_and(Stellar_Mass > (9.0),DustRate_GROW>0.0)
    DustRate_GROW_z[loop]= np.log10(np.sum(DustRate_GROW[condition==1])/volume)
    condition = np.logical_and(Stellar_Mass > (9.0),DustRate_DEST>0.0)
    DustRate_DEST_z[loop]= np.log10(np.sum(DustRate_DEST[condition==1])/volume)
    condition = np.logical_and(Stellar_Mass > (9.0),DustRate_ALL>0.0)
    DustRate_ALL_z[loop] = np.log10(np.sum(DustRate_ALL[condition==1]-DustRate_DEST[condition==1])/volume)
    condition = np.logical_and(Stellar_Mass > (9.0),DustRate_ALL>0.0)
    SFR_z[loop] = np.log10(np.sum(SFR[condition==1])/volume)
    


    ##print DustRate_z[loop],DustRate_AGB_z[loop],DustRate_SNII_z[loop],DustRate_SNIA_z[loop],DustRate_GROW_z[loop],DustRate_DEST_z[loop]



    if loop==0:
        #fig = plt.figure(figsize=(9,9))
        fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(9,9))
        fig.subplots_adjust(hspace=0)
        fig.subplots_adjust(wspace=0)
    plt.subplot(3,3,loop+1)
    plt.xlim([6,11.98])
    plt.ylim([-10,2.47])
    if(loop==1 or loop==2 or loop==4 or loop==5):
        plt.tick_params(axis='both', which='major', labelsize=12,width=2,length=6, labelleft ='off', labelbottom='off', direction='inout')
        plt.tick_params(axis='both', which='minor', labelsize=12,width=2,length=6, labelleft ='off', labelbottom='off', direction='inout')

    if(loop==0 or loop==3):
        plt.tick_params(axis='both', which='major', labelsize=12,width=2,length=6, labelleft ='on', labelbottom='off', direction='inout')
        plt.tick_params(axis='both', which='minor', labelsize=12,width=2,length=6, labelleft ='on', labelbottom='off', direction='inout')

    if(loop==7 or loop==8):
        plt.tick_params(axis='both', which='major', labelsize=12,width=2,length=6, labelleft ='off', labelbottom='on', direction='inout')
        plt.tick_params(axis='both', which='minor', labelsize=12,width=2,length=6, labelleft ='off', labelbottom='on', direction='inout')
    
    if(loop==6):
        plt.tick_params(axis='both', which='major', labelsize=12,width=2,length=6, labelleft ='on', labelbottom='on', direction='inout')
        plt.tick_params(axis='both', which='minor', labelsize=12,width=2,length=6, labelleft ='on', labelbottom='on', direction='inout')
    
    if(loop==0):
        plt.ylim([-10,2.5])
    if(loop==8):
        plt.xlim([6,12])


    if loop==7:
        plt.xlabel(r'log$_{10}$(M$_{*}$/M$_{\odot}$)', fontsize=18)
    if loop==3:
        plt.ylabel(r'Dust production rate M$_{\odot}$yr$^{-1}$', fontsize=18)
    #plt.tick_params(axis='both', which='major', labelsize=12,width=2)
    #plt.tick_params(axis='both', which='minor', labelsize=12,width=2)
    #plt.text(7,1,"N = "+str(sum(count))+"\nz = "+str(loop))
    plt.text(9,-9,"z = "+str(loop),fontsize=16)
    
    plt.errorbar(SM_bins_AGB ,Dust_med_AGB  ,yerr=(Dust_mederr_AGB) ,color='b',label='AGB',linewidth=2)
    plt.errorbar(SM_bins_SNII,Dust_med_SNII ,yerr=(Dust_mederr_SNII),color='r',label='SNII',linewidth=2)
    plt.errorbar(SM_bins_SNIA,Dust_med_SNIA ,yerr=(Dust_mederr_SNIA),color='y',label='SNIA',linewidth=2)
    plt.errorbar(SM_bins_GROW,Dust_med_GROW ,yerr=(Dust_mederr_GROW),color='g',label='GG',linewidth=2)
    plt.errorbar(SM_bins_DEST,Dust_med_DEST ,yerr=(Dust_mederr_DEST),color='k',label='DEST',linewidth=2)
    #plt.errorbar(SM_bins_ALL ,Dust_med_ALL  ,yerr=(Dust_mederr_ALL),color='orange',label='ALL')

axes = fig.get_axes()
for ax in axes:
    [i.set_linewidth(2.1) for i in ax.spines.itervalues()]


pylab.savefig('./graphs/stellar_rate_ALL.png', bbox_inches=0)
plt.close()
    







fig = plt.figure(figsize=(7,7))
plt.xlabel(r'redshift', fontsize=18,labelpad=10)
plt.ylabel(r'log$_{10}$(Cosmic dust rate Msol/yr/Mpc$^3$)', fontsize=18,labelpad=0)
plt.xlim(-1,10)

plt.tick_params(axis='both', which='major', labelsize=12,width=2)
plt.tick_params(axis='both', which='minor', labelsize=12,width=2)

plt.plot(DustRate_z,DustRate_AGB_z , color='b',label='AGB',linewidth = 2)
plt.plot(DustRate_z,DustRate_SNII_z, color='r',label='SNII',linewidth = 2)
plt.plot(DustRate_z,DustRate_SNIA_z, color='y',label='SNIA',linewidth = 2)
plt.plot(DustRate_z,DustRate_GROW_z, color='g',label='GG',linewidth = 2)
plt.plot(DustRate_z,DustRate_DEST_z, color='k',label='DEST',linewidth = 2)
plt.plot(DustRate_z,DustRate_ALL_z,  color='orange',label='ALL',linewidth = 2)
#plt.plot(DustRate_z,DustRate_GROW_z+DustRate_AGB_z+DustRate_SNIA_z+DustRate_SNII_z-DustRate_DEST_z, color='cyan', linestyle='--', linewidth=2,label='Net')
plt.plot(DustRate_z,SFR_z, color='orange', linestyle='--', linewidth=2,label='SFR')
plt.legend(loc='lower right')

axes = fig.get_axes()
for ax in axes:
    [i.set_linewidth(2.1) for i in ax.spines.itervalues()]


pylab.savefig('./graphs/rate_redshift.png', bbox_inches=0)
plt.close()













