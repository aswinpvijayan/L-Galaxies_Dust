#=========================================================================
#
#  Script to read in L-galaxies snapshot data
#
#  To force a re-read of the data do Gal=None
#
#-------------------------------------------------------------------------

# Imports

import sys

datadir = '../../output/'

sys.path.insert(0,datadir)

# Template structure for L-Galaxies data
import snap_template   # structure temple for data
import read_lgal       # function to read in data

#-------------------------------------------------------------------------

# Parameters

# Snaplist file
snaplist_file = '../MRPlancksnaplist.txt'

# Define what snapshot you want



for i in range(0,10):
	if i==0:
		snapshot=58
		file_prefix = "SA_z0.00"
		output_file = "../data/lgal_z0.pkl"
	elif i==1:
		snapshot=38
		file_prefix = "SA_z1.04"
		output_file = "../data/lgal_z1.pkl"
	elif i==2:
		snapshot=30
		file_prefix = "SA_z2.07"
		output_file = "../data/lgal_z2.pkl"
	elif i==3:
		snapshot=25
		file_prefix = "SA_z3.11"
		output_file = "../data/lgal_z3.pkl"
	elif i==4:
		snapshot=22
		file_prefix = "SA_z3.95"
		output_file = "../data/lgal_z4.pkl"
	elif i==5:
		snapshot=19
		file_prefix = "SA_z5.03"
		output_file = "../data/lgal_z5.pkl"
	elif i==6:
		snapshot=17
		file_prefix = "SA_z5.92"
		output_file = "../data/lgal_z6.pkl"
	elif i==7:
		snapshot=15
		file_prefix = "SA_z6.97"
		output_file = "../data/lgal_z7.pkl"
	elif i==8:
		snapshot=13
		file_prefix = "SA_z8.22"
		output_file = "../data/lgal_z8.pkl"
	elif i==9:
		snapshot=12
		file_prefix = "SA_z8.93"
		output_file = "../data/lgal_z9.pkl"


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
	props['DustMassISM'] = True
	props['DustMassCGM'] = True

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
	fout = open(output_file,'wb')
	cPickle.dump(gals,fout,cPickle.HIGHEST_PROTOCOL)
	fout.close()
