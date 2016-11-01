# Script to read in pickled file of galaxy properties.
# Requires pickleFile to hold name of desired input file.
# Returns galaxy properties in gals.

import cPickle

fin = open('../data/lgal_output.pkl','rb')
gals=cPickle.load(fin)
fin.close()

print gals['Type']
print gals['DustMass']

print gals['Type'][0]
print gals['DustMass'][0][1]

# print len(gals['Type'][0])
# print len(gals['DustMass'][0][0])
# print len(gals['DustMass'])/24

# for i in range(0,
