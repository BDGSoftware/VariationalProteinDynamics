
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
print ("vp_basisVectors.py")
print ("Calculates residue-centered basis vectors for three-state microstate trajectories")
print ("---------------------------------------------------") 
print "START: "+str(datetime.now())
print ("---------------------------------------------------")

#------------------------------------------------
# Command line arguments  
#------------------------------------------------
parser = argparse.ArgumentParser(description='Calculates residue-centered basis vectors for three-state microstate trajectories')
parser.add_argument('-inDir','--inDir', help='Directory with microstate-trajectories (npy)',required=True)
parser.add_argument('-outDir','--outDir', help='Directory to which the basis vectors are written (npy)',required=True)
parser.add_argument('-sequence','--sequence', help='File with peptide sequence and included residues',required=True)
parser.add_argument('-basisSet','--basisSet', help='Directory with residue centered basis vectors to be orthonormalized',required=False, default=False)
args = parser.parse_args()

print ("Using the following command line arguments")
print ("-inDir (Directory with microstate-trajectories (npy)):")
print ("%s" % args.inDir) 
print ("-outDir (Directory to which the basis vectors are written (npy):")
print ("%s" % args.outDir)
print ("-sequence (File with peptide sequence and included residues):")
print ("%s" % args.sequence)
print ("-basisSet (Directory with residue centered basis vectors to be orthonormalized):")
add_remark=''
if args.basisSet==False:
	add_remark='. Internal standard basis functions will be used.'
print ("%s" % args.basisSet+ add_remark)
print ("---------------------------------------------------")

#------------------------------------------------
# Functions  
#------------------------------------------------
def scp(v,w,pi):# weighted scalar product
	return np.multiply(v,pi).dot(w)

def norm(v,pi):# weighted norm
	return np.sqrt(scp(v,v,pi))

def gram_schmidt(basis,pi):# gram-schmidt method
	on_basis=np.zeros(basis.shape)
	for i,v in enumerate(basis.transpose()):
		v_sc=np.array(map(lambda x: scp(x,v,pi), list(on_basis.transpose())))
		e=v-np.dot(on_basis,v_sc)
		e=e/norm(e,pi)
		on_basis[:,i]=e
	return on_basis

#------------------------------------------------
# Internal basis sets  
#------------------------------------------------
RBV=np.array([[1,1,1],[-1,-1,1],[-1,1,0]]).transpose()
RBVp=np.array([[1,1],[-1,1]]).transpose()
RBVg=np.array([[1,1,1,1,1,1],[1,-1,1,-1,1,0],[-1,-1,1,1,1,-1]]).transpose()
RBVxp=np.array([[1,1,1],[-1,1,1],[-1,1,-1]]).transpose()
RBVgp=np.array([[1,1,1,1,1,1],[1,-1,1,-1,0,1],[1,-1,-1,-1,0,1]]).transpose()
#------------------------------------------------
# Read sequence file  
#------------------------------------------------
sequence=[]
# open file
with open(args.sequence, 'r') as f: 
    #loop over lines
    for line in f:
        # check wether the line is a comment        
        if not line.startswith("#"):
            sequence.append(line.strip().split("\t"))

# get list of included residues
# extract the second column of sequence and convert to int
resIncluded=[int(i[1]) for i in sequence] 
# convert to numpy array
resIncluded = np.array(resIncluded)

# get the total number of included residues
nResIncluded = resIncluded.sum()
# extract the actual sequence from the array sequence
sequence=np.array(zip(*sequence)[0])

#------------------------------------------------
# Compute basis vectors
#------------------------------------------------
print "Reading microstate trajectories"

for i, res in enumerate(sequence):
    if resIncluded[i]==1:

	# identify rbv basis
	prep=i<len(sequence)-1 and sequence[i+1]=='P' and res!='P'
	if args.basisSet==False:
		if res=='P':
			basis=RBVp	
		elif res=='G':
			if prep:
				basis=RBVgp
			else:
				basis=RBVg	
		else:
			if prep:
				basis=RBVxp
			else:
				basis=RBV
	else:
		if prep:
			basis=np.loadtxt(args.basisSet+'/Ac_'+res+'P_NHMe/rEV_norm.dat')
		else:
			basis=np.loadtxt(args.basisSet+'/Ac_'+res+'_NHMe/rEV_norm.dat')	
	n_ms=basis.shape[0]

        # construct paths to the microstate trajectory files
	microstateTrajFiles=[]
	for root,dirs,files in os.walk(args.inDir):
		for file in files:
			if file==str(i+1)+res+'.npy':
				fname = os.path.join(root,file)
				microstateTrajFiles.append(np.load(fname))
      	
	# concatenate trajectories
	mtraj=np.hstack(microstateTrajFiles)
	
	# compute equillibrium distribution
	pi=np.bincount(mtraj)/float(mtraj.shape[0])
	p_pop=np.where(pi!=0)[0].shape[0]*100/float(pi.shape[0])	
	
	# print useful information to the terminal
	print res+ ':\t Using a '+str(n_ms) +' state discretization.'
	print '\t Found '+str(len(microstateTrajFiles))+' microstate trajectories.'
	print '\t '+str(round(p_pop,1))+' % of microstates populated.'

	# Check validity of data dimensions and give warning if necessary
	if pi.shape[0]!=n_ms:
		
		diff=basis.shape[0]-pi.shape[0]
		print 'WARNING: The last '+str(diff)+' microstates are not populated. This could indicate mismatched input data.'
		opt = bool(input('Please choose:\n 0\t Abort program\n 1\t Ignore and continue (will append zeros)\n'))
	
		# Fix data dimensions
		if opt:
			pi=np.hstack((pi,np.zeros(diff)))
		else:
			sys.exit(1)

	# calculate ON basis
	try:
		basis=gram_schmidt(basis,pi)
	except ValueError:
		print 'ERROR: Microstate discretization and basis vector dimension do not match. Aborting.'
		sys.exit(1)

	# write basis into file
	np.savetxt(args.outDir+'/'+str(i+1)+res+'.dat',basis)
	print '\t Saving orthonormal basis to "'+args.outDir+'/'+str(i+1)+res+'.dat".'







