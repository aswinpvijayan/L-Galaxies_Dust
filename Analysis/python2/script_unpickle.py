# Script to read in pickled file of galaxy properties.
# Requires pickleFile to hold name of desired input file.
# Returns galaxy properties in gals.

import cPickle

fin = open(pickleFile, 'rb')
gals=cPickle.load(fin)
fin.close()
