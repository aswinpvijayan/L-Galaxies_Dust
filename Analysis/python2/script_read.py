#=========================================================================
#
#  Script to read in L-galaxies snapshot data
#
#  To force a re-read of the data do Gal=None
#
#-------------------------------------------------------------------------

# Imports

import sys

# Path to data
#datadir = '/mnt/lustre/scratch/virgo/SAM_output/Hen14_sfh2'
#datadir = '/mnt/lustre/scratch/bf77/L-Galaxies/Hen14_MR/output'
#datadir = '/lustre/scratch/astro/bf77/L-Galaxies/Hen14_MR/test_run'
#datadir = '/lustre/scratch/astro/bf77/L-Galaxies/Hen14_MR/output_ref'
#datadir = '/mnt/lustre/scratch/bf77/L-Galaxies/Hen14_MR/output_XS'
#datadir = '/mnt/lustre/scratch/bf77/L-Galaxies/Hen14_MR/output_cooling_new'

datadir = '../../output/'
#datadir = '/lustre/scratch/astro/bf77/HANNAH_test/'
#datadir = '/lustre/scratch/astro/bf77/L-Galaxies/Hen14_MR/output_newcoolx8'
#datadir = '/lustre/scratch/astro/bf77/L-Galaxies/Hen14_MR/output_newcoolx10infmod' #test_run'

#datadir = '/mnt/lustre/scratch/bf77/output_test'
#datadir = '/mnt/lustre/scratch/bf77/L-Galaxies/Hen14_MR/output_test_sample'
sys.path.insert(0,datadir)

# Template structure for L-Galaxies data
import snap_template   # structure temple for data
import read_lgal       # function to read in data

#-------------------------------------------------------------------------

# Parameters

# Snaplist file
snaplist_file = '../MRPlancksnaplist.txt'

# Define what snapshot you want
snapshot = 58 #58 => z = 0

# Define which files you want to read in
firstfile = 5
lastfile = 5 #511

# Define what properties you want to read in
props = snap_template.properties_used

props['Type'] = True
props['ColdGas'] = True
props['StellarMass'] = True
props['BulgeMass'] = True
props['DiskMass'] = True
props['HotGas'] = True
props['ICM'] = True
props['MetalsColdGas'] = True
props['MetalsBulgeMass'] = True
props['MetalsDiskMass'] = True
props['MetalsHotGas'] = True
props['MetalsEjectedMass'] = True
props['MetalsICM'] = True
props['Sfr'] = True
props['SfrBulge'] = True
props['DiskMass_elements'] = True
props['BulgeMass_elements'] = True
props['ColdGas_elements'] = True
props['DustMass'] = True
# 

#-------------------------------------------------------------------------

# Working body of the program

# Read in redshift of snapshot and create file prefix
#f = open(snaplist_file)
#lines = f.readlines()
#f.close()
#for this_line in lines:
#    words = this_line.split()
    #print words[0],words[2]
#    if words[0]==str(snapshot):
#        file_prefix = "SA_z"+words[2]
#        print "z = ", words[2]
file_prefix = "SA_z0.00"
# Read in galaxy output
(nTrees,nHalos,nTreeHalos,gals) = \
    read_lgal.read_snap(datadir,file_prefix,firstfile,lastfile,\
                            props,snap_template.struct_dtype)


import cPickle

#fout = open('/mnt/lustre/scratch/bf77/Dave_lum.pkl', 'wb')
#fout = open('/mnt/lustre/scratch/bf77/Lgal_coolXS_less.pkl', 'wb')
#fout = open('/lustre/scratch/astro/bf77/Lgal_runx10mcmc_z0_part1.pkl', 'wb')
#fout = open('/lustre/scratch/astro/bf77/Lgal_runx10infmod_z0_part1.pkl', 'wb')
#fout = open('/lustre/scratch/astro/bf77/try3.pkl', 'wb') #Lgal_runx10AGNXS01_z0_part1.pkl', 'wb')
#fout = open('/lustre/scratch/astro/bf77/Lgal_run_X10_NewReinc_NoQuasar_Infall0_less.pkl', 'wb')
fout = open('../data/lgal_output.pkl','wb')
cPickle.dump(gals,fout,cPickle.HIGHEST_PROTOCOL)
fout.close()
