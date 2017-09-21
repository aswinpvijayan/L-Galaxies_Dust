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


mag_clay15_int_z4, phi_clay15_int_z4 = np.loadtxt('/Users/scottclay/Desktop/apollo/Clay15_z4.txt',unpack=True,comments='#')
mag_clay15_int_z5, phi_clay15_int_z5 = np.loadtxt('/Users/scottclay/Desktop/apollo/Clay15_z5.txt',unpack=True,comments='#')
mag_clay15_int_z6, phi_clay15_int_z6 = np.loadtxt('/Users/scottclay/Desktop/apollo/Clay15_z6.txt',unpack=True,comments='#')
mag_clay15_int_z7, phi_clay15_int_z7 = np.loadtxt('/Users/scottclay/Desktop/apollo/Clay15_z7.txt',unpack=True,comments='#')

mag_clay15_obs_z4, phi_clay15_obs_z4 = np.loadtxt('/Users/scottclay/Desktop/apollo/Clay15_z4_obs.txt',unpack=True,comments='#')
mag_clay15_obs_z5, phi_clay15_obs_z5 = np.loadtxt('/Users/scottclay/Desktop/apollo/Clay15_z5_obs.txt',unpack=True,comments='#')
mag_clay15_obs_z6, phi_clay15_obs_z6 = np.loadtxt('/Users/scottclay/Desktop/apollo/Clay15_z6_obs.txt',unpack=True,comments='#')
mag_clay15_obs_z7, phi_clay15_obs_z7 = np.loadtxt('/Users/scottclay/Desktop/apollo/Clay15_z7_obs.txt',unpack=True,comments='#')


for loop in range(4,5):
#for loop in range(0,1):

    fin = open('../data/MR/lgal_z'+str(loop)+'.pkl','rb')
    gals_MR=cPickle.load(fin)
    fin.close()

    fin = open('../data/MR_UV/lgal_z'+str(loop)+'.pkl','rb')
    gals_UV=cPickle.load(fin)
    fin.close()
    ##############
    

    Stellar_Mass_MR = np.zeros(len(gals_MR['Type']))
    Dust_Mass_MR = np.zeros(len(gals_MR['Type']))
    MAG_MR = np.zeros(len(gals_MR['Type']))
    MAGDUST_MR = np.zeros(len(gals_MR['Type']))
    Dust_Mass = np.zeros(len(gals_MR['Type']))
    Metal_Mass = np.zeros(len(gals_MR['Type']))
    ColdGas = np.zeros(len(gals_MR['Type']))
    GasDiskRadius = np.zeros(len(gals_MR['Type']))
    
    
    Stellar_Mass_MR = gals_MR['StellarMass']*1.0E10/0.673
    Dust_Mass_MR = np.sum(gals_MR['Dust_elements'], axis=1)
    ColdGas = gals_MR['ColdGas']*1.0E10/0.673
    GasDiskRadius = gals_MR['GasDiskRadius']

    MAG_MR = gals_MR['Mag'][:,0]
    MAGDUST_MR = gals_MR['MagDust'][:,0]
    
    MAG_UV = gals_UV['Mag'][:,0]
    MAGDUST_UV = gals_UV['MagDust'][:,0]
    
    

#     Dust_Mass = np.sum(gals_MR['Dust_elements'], axis=1)
#     Metal_Mass= np.sum(gals_MR['ColdGas_elements'][:,2:11])
#     
#     DTM = Dust_Mass / Metal_Mass
#     
#     print (np.log10(DTM))
#     
# 
#     N_H = (ColdGas/(1.0E10/0.673)) / (3.14159 * (GasDiskRadius * 0.94)**2 * 1.4)
#     N_H /= 3252.37#;    // 3252.37 = 10^(3.5122) ... ha ha ! 
# 
#     
#     A_V = 0.45 * (DTM/0.737) * (((Dust_Mass + Metal_Mass)/ColdGas)/0.0137) * (N_H)
#     
#     A_1500 = A_V * 3.6 # 540nm/1500A = 3.6
#     
#     MAG_NEW = MAG_MR + A_1500
#     
#     print "--------------"
#     
#     print A_V
#     print (DTM/0.737)
#     print (((Dust_Mass + Metal_Mass)/ColdGas)/0.0137)
#     print (N_H)
# 
# 
#     print "--------------"
    
    #for i in range(0,len(gals_MR['Type'])):
    #    print MAG_MR[i], A_1500[i], MAG_NEW[i], (DTM[i]/0.737),N_H[i]
    
    
    #---------------------UVLF function
    #volume_MR = (480.279/0.673)**3.0
    
    condition = np.logical_and(MAG_MR>-90.0, MAG_MR<-12.0)
    UV_MR = MAG_MR[condition==1] - 5*np.log10(0.673)
    hist_MR,bin_edges_MR = np.histogram(UV_MR,bins=50)
    bin_centers_MR = (bin_edges_MR[:-1] + bin_edges_MR[1:]) / 2
    #volume_MR = ((480.279/512)/0.673)**3.0
    volume_MR = ((480.279)**3) * 1/512
    binsize_MR = (bin_centers_MR[-1] - bin_centers_MR[0])/len(bin_centers_MR)
    
    
    condition = np.logical_and(MAGDUST_MR>-90.0, MAGDUST_MR<-12.0)
    UV_MR_DUST = MAGDUST_MR[condition==1]- 5*np.log10(0.673)
    hist_MR_DUST,bin_edges_MR_DUST = np.histogram(UV_MR_DUST,bins=50)
    bin_centers_MR_DUST = (bin_edges_MR_DUST[:-1] + bin_edges_MR_DUST[1:]) / 2
    #volume_MR = ((480.279/512)/0.673)**3.0
    volume_MR = ((480.279)**3) * 1/512
    binsize_MR_DUST = (bin_centers_MR_DUST[-1] - bin_centers_MR_DUST[0])/len(bin_centers_MR_DUST)
    
    condition = np.logical_and(MAG_UV>-90.0, MAG_UV<-12.0)
    UV_UV = MAG_UV[condition==1] - 5*np.log10(0.673)
    hist_UV,bin_edges_UV = np.histogram(UV_UV,bins=50)
    bin_centers_UV = (bin_edges_UV[:-1] + bin_edges_UV[1:]) / 2
    volume_UV = ((480.279)**3) * 1/512
    binsize_UV = (bin_centers_UV[-1] - bin_centers_UV[0])/len(bin_centers_UV)
    
    condition = np.logical_and(MAGDUST_UV>-90.0, MAGDUST_UV<-12.0)
    UV_UV_DUST = MAGDUST_UV[condition==1]- 5*np.log10(0.673)
    hist_UV_DUST,bin_edges_UV_DUST = np.histogram(UV_UV_DUST,bins=50)
    bin_centers_UV_DUST = (bin_edges_UV_DUST[:-1] + bin_edges_UV_DUST[1:]) / 2
    volume_UV = ((480.279)**3) * 1/512
    binsize_UV_DUST = (bin_centers_UV_DUST[-1] - bin_centers_UV_DUST[0])/len(bin_centers_UV_DUST)
    
    
    #plt.xlim([-25,-15])
    #plt.ylim([0,9.98])

    plt.plot(mag_clay15_int_z4,np.log10(phi_clay15_int_z4), color='k', label = 'Clay15 Intrinsic', linewidth = 2)
    plt.plot(mag_clay15_obs_z4,np.log10(phi_clay15_obs_z4), color='k',linestyle='--', label = 'Clay15 Observed', linewidth = 2)
    
    plt.plot(bin_centers_MR,np.log10(hist_MR/(volume_MR*binsize_MR)),linestyle='-',color='r',label='Dust int',linewidth=2)
    plt.plot(bin_centers_MR_DUST,np.log10(hist_MR_DUST/(volume_MR*binsize_MR_DUST)),linestyle='--',color='red',label='Dust obs',linewidth=2)
    #plt.plot(bin_centers_NEW,np.log10(hist_NEW/(volume_NEW*binsize_NEW)),linestyle='-',color='blue',label='New Median',linewidth=2)

    plt.plot(bin_centers_UV,np.log10(hist_UV/(volume_UV*binsize_UV)),linestyle=':',color='b',label='UV Int',linewidth=2)
    plt.plot(bin_centers_UV_DUST,np.log10(hist_UV_DUST/(volume_UV*binsize_UV_DUST)),linestyle=':',color='b',label='UV obs',linewidth=2)

    plt.legend()
    pylab.savefig('./graphs/UVLF_z4.png', bbox_inches=0)
    plt.close()
    
#     plt.xlim([6,12])
#     #plt.ylim([-8,2])
#     
#      volume = (480.279/0.673)**3.0
# #     volume = (0.938044921/0.673)**3.0
#     binsize = 1.35135751
#     
#     plt.plot(bin_centers,hist/(volume*binsize),color='k',label='Dust All')
#     
#     plt.xlabel(r'log$_{10}$(Mdust/M$_{\odot}$)', fontsize=14,labelpad=10)
#     #plt.ylabel(r'log$_{10}$(Mdust/Mcoldgas$)', fontsize=14,labelpad=0)
#     plt.tick_params(axis='both', which='major', labelsize=10)
#     plt.tick_params(axis='both', which='minor', labelsize=8)
#     #plt.legend(loc='lower right')
#     #plt.text(2,-7,"N = "+str(sum(count)))
#     #plt.text(2,-5,"z = "+str(loop)+" :All dust")
#     
#     pylab.savefig('./graphs/dustmass_function_z'+str(loop)+'.png', bbox_inches=0)
#     plt.close()
# 


















