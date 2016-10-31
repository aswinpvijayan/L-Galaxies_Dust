# script to call other scripts to read in, mass trim and pickle the data

#model='Hen14'
#model='Hen14_spin'
#model='Guo10_spin'
model='Guo10'
#model='Hen14_MRII'

for redshift in (0,0.15,0.25,0.5,1,1.5,2,3,4,5,6):
#for redshift in (1.5,2,3,4,5,6):
    print('Starting redshift ',redshift)
    execfile('script_read.py')
    execfile('script_masstrim.py')
    execfile('script_pickle.py')
