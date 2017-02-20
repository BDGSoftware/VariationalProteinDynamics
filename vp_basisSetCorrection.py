
#------------------------------------------------
#  Imports
#------------------------------------------------
import argparse
import numpy as np
import os
import sys
from datetime import datetime

# output to screen  
print ("---------------------------------------------------")
print ("vp_3basisSetCorrection.py")
print ("Corrects the basis set based on deviations from orthonormality and corrects correlation matrices accordingly.")
print ("---------------------------------------------------") 
print "START: "+str(datetime.now())
print ("---------------------------------------------------")

#------------------------------------------------
# Command line arguments  
#------------------------------------------------
parser = argparse.ArgumentParser(description='Corrects the basis set based on deviations from orthonormality and corrects correlation matrices accordingly.')
parser.add_argument('-inDir','--inDir', help='Directory with correlation matrices',required=True)
parser.add_argument('-outDir','--outDir', help='Directory to which the corrected correlation matrices are written',required=True)
parser.add_argument('-PBVFile','--PBVFile', help='Basis set file',required=True)
parser.add_argument('-PBVFileOut','--PBVFileOut', help='Path to the corrected basis set file to be written',required=True)

args = parser.parse_args()

print ("Using the following command line arguments")
print ("-inDir (Directory with correlation matrices):")
print ("%s" % args.inDir) 
print ("-outDir (Directory to which the corrected correlation matrices are written:")
print ("%s" % args.outDir)
print ("-PBVFile (Basis set file):")
print ("%s" % args.PBVFile)
print ("-PBVFileOut (Path to the corrected basis set file to be written):")
print ("%s" % args.PBVFileOut)
print ("---------------------------------------------------")

#------------------------------------------------
# Correct basis set  
#------------------------------------------------
# load overlap diagonal 
SDiag=np.diag(np.load(args.inDir+'/S.npy'))

# get indices to be excluded
i=np.nonzero((SDiag<0.1)+(SDiag>10))[0]

# print message
print 'The following basis functions will be excluded: '
print i

#------------------------------------------------
# Create new basis set file 
#------------------------------------------------
# load old basis set 
PBVFile=np.loadtxt(args.PBVFile)

# create new basis set
PBVFileNew=np.delete(PBVFile,i,axis=0)

# save new basis set
np.savetxt(args.PBVFileOut, PBVFileNew, delimiter='\t', fmt='%i')

#------------------------------------------------
# Correct correlation matrices
#------------------------------------------------
# identify matrices to be corrected

matrices=[]
matricesIdentifier=[]
for root, dirs, files in os.walk(args.inDir):
	for file in files:
		if file.endswith('.npy'):
			fname =os.path.join(root,file)
			matrices.append(np.load(fname))
			matricesIdentifier.append(file)

# loop over matrices
for k,m in enumerate(matrices):
	print 'correcting '+ matricesIdentifier[k]
	
	# delete rows and columns
	mC=np.delete(m,i,0)
	mC=np.delete(mC,i,1)
	
	# save new matrix
	np.save(args.outDir+'/'+matricesIdentifier[k],mC)





