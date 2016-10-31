#=========================================================================
#
#  Script to read in L-galaxies snapshot data
#
#  To force a re-read of the data do gals=None
#
#-------------------------------------------------------------------------

# Imports

import sys

# Template structure for L-Galaxies data
import read_lgal       # function to read in data

#-------------------------------------------------------------------------

# Parameters

# Decide what data set you want
# Allow for calling with model set from metascript
try:
    model
except:
    #model='Hen14'
    #model='Hen14_spin'
    model='Guo10'
    #model='Guo10_spin'
    #model='Hen14_MRII'

# Decide what redshift you want
# Allow for calling with redshift set from metascript
try:
    redshift
except:
    redshift=0.15

# Define which files you want to read in
# 0-511 gives all the files (but takes a long time)
firstfile = 0
lastfile = 511

# Path to input data
if model=='Hen14': idatadir = '/lustre/scratch/astro/virgo/SAM_output/Hen14_sfh2'
if model=='Hen14_spin': idatadir = '/lustre/scratch/astro/virgo/SAM_output/Hen14_spin'
if model=='Guo10': idatadir = '/lustre/scratch/astro/virgo/SAM_output/Guo10_sfh2'
if model=='Guo10_spin': idatadir = '/lustre/scratch/astro/virgo/SAM_output/Guo10_spin'
if model=='Hen14_MRII': idatadir = '/lustre/scratch/astro/virgo/SAM_output/Hen14_MRII'

# Path for output data
if model=='Hen14': odatadir = 'data/Hen14_sfh2'
if model=='Hen14_spin': odatadir = 'data/Hen14_spin'
if model=='Guo10': odatadir = 'data/Guo10_sfh2'
if model=='Guo10_spin': odatadir = 'data/Guo10_spin'
if model=='Hen14_MRII': odatadir = 'data/Hen14_MRII'

# import template
sys.path.insert(0,idatadir)
import snap_template

# Define what properties you want to read in
props = snap_template.properties_used
props['Type'] = True
props['Rvir'] = True
props['Mvir'] = True
props['CentralMvir'] = True
props['DistanceToCentralGal'] = True
props['Pos'] = True
props['Vel'] = True
props['Vmax'] = True
props['InfallVmax'] = True
props['InfallVmaxPeak'] = True
if model=='Hen14_spin':
    props['BulgeSpin'] = True
    props['BulgeSpinMax'] = True
    props['DiskSpin'] = True
    props['DiskSpinMax'] = True
    props['ColdGasSpin'] = True
props['BulgeMass'] = True
props['DiskMass'] = True
props['StellarMass'] = True
props['ColdGas'] = True
props['HotGas'] = True
props['EjectedMass'] = True
props['MetalsDiskMass'] = True
props['MetalsStellarMass'] = True
props['MetalsBulgeMass'] = True
props['MetalsColdGas'] = True
props['MetalsHotGas'] = True
props['MetalsEjectedMass'] = True
props['BulgeSize'] = True
if model=='Hen14_spin':
    props['DiskRadius'] = True
    props['ColdGasRadius'] = True
props['Sfr'] = True
props['BlackHoleMass'] = True
props['BlackHoleGas'] = True
props['QuasarAccretionRate'] = True
props['RadioAccretionRate'] = True
if model=='Hen14_spin' or model=='Guo10_spin':
    props['ObsMag'] = True
else:
    props['Mag'] = True
props['sfh_ibin'] = True
props['sfh_numbins'] = True
props['sfh_DiskMass'] = True
props['sfh_BulgeMass'] = True
props['sfh_MetalsDiskMass'] = True
props['sfh_MetalsBulgeMass'] = True


#-------------------------------------------------------------------------

# Working body of the program

# Matching between redshift and snapshot provided by dictionary
# In future create this dictionary from appropriate files
snapz_Hen14={0:58,0.15:53,0.25:50,0.5:45,1:38,1.5:34,2:30,2.5:28,3:25,4:22,5:19,6:17}
snapz_Guo10={0:63,0.15:57,0.25:54,0.5:48,0.75:52,1:41,1.5:36,2:32,3:27,4:24,5:21,6:18}
snapz_Hen14_MRII={0:62,0.25:54,0.5:49,1:42,1.5:38,2:34,3:29,4:26,5:23,6:21}
if model=='Hen14': snapshot = snapz_Hen14[redshift]
if model=='Hen14_spin': snapshot = snapz_Hen14[redshift]
if model=='Guo10': snapshot = snapz_Guo10[redshift]
if model=='Guo10_spin': snapshot = snapz_Guo10[redshift]
if model=='Hen14_MRII': snapshot = snapz_Hen14_MRII[redshift]

# Snaplist file
if model=='Hen14': snaplist_file = idatadir+'/MRPlancksnaplist.txt'
if model=='Hen14_spin': snaplist_file = idatadir+'/MRPlancksnaplist.txt'
if model=='Guo10': snaplist_file = idatadir+'/MRW1snaplist.txt'
if model=='Guo10_spin': snaplist_file = idatadir+'/MRW1snaplist.txt'
if model=='Hen14_MRII': snaplist_file = idatadir+'/MRIIPlancksnaplist.txt'

# Read in redshift of snapshot and create file prefix
f = open(snaplist_file)
lines = f.readlines()
f.close()
for this_line in lines:
    words = this_line.split()
    #print words[0],words[2]
    if words[0]==str(snapshot):
        file_prefix = "SA_z"+words[2]

# Read in galaxy output
(nTrees,nHalos,nTreeHalos,gals) = \
    read_lgal.read_snap(idatadir,file_prefix,firstfile,lastfile,\
                            props,snap_template.struct_dtype)
