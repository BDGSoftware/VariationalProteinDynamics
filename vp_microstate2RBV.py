# -*- coding: utf-8 -*-
"""
Created on Fri May 22 09:19:56 2015

@author: bettina
"""

#------------------------------------------------
#  Imports
#------------------------------------------------
import argparse
import numpy as np
import os
import sys
from datetime import datetime

#--------------------------------
__author__ = 'bettina'

# output to screen  
print ("---------------------------------------------------")
print ("vp_microstate2PBV.py")
print ("Projects microstate-trajectories onto the peptide-based vectors (PBVs)")
print ("B. Keller (bettina.keller@fu-berlin.de), May. 2015")
print ("---------------------------------------------------") 
print "START: "+str(datetime.now())
print ("---------------------------------------------------")

#------------------------------------------------
# Command line arguments  
#------------------------------------------------
parser = argparse.ArgumentParser(description='Projects microstate-trajectories onto the peptide-based vectors (PBVs)')
parser.add_argument('-inDir','--inDir', help='Directory with microstate-trajectories (npy)',required=True)
parser.add_argument('-outDir','--outDir', help='Directory to which the PBV trajectories are written (npy)',required=True)
parser.add_argument('-sequence','--sequence', help='File with peptide/protein sequence and included residues',required=True)
parser.add_argument('-basisSet', '--basisSet', default=os.path.dirname(sys.argv[0])+'/BasisVectors', help='Directory with basis set. Default: ./BasisVectors')
args = parser.parse_args()


print ("Using the following command line arguments")
print ("-inDir (Directory with microstate-trajectories (npy)):")
print ("%s" % args.inDir) 
print ("-outDir (Directory with microstate-trajectories (npy):")
print ("%s" % args.outDir)
print ("-sequence (File with peptide/protein sequence and included residues):")
print ("%s" % args.sequence)
print ("-basisSet (Directory with basis set. Default: ./BasisVectors):")
print ("%s" % args.basisSet)
print ("---------------------------------------------------")

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
# Read residue-based vectors (RBVs)  
#------------------------------------------------
print "Reading residue-based vectors"
# List of RBVs
RBV=[]
# loop over residues
# the RBVs of each residue are added to the list
# if a given residue type occurs several times, the corresponding RBVs
# are read and added to the list every time the residue is marked as included
for i, res in enumerate(sequence):
	if resIncluded[i]==1:
		if os.path.exists(args.basisSet+'/'+str(i+1)+res+'.dat'):
			RBVFile=args.basisSet+'/'+str(i+1)+res+'.dat'
		elif i+1<len(sequence) and sequence[i+1]=='P' and sequence[i]!='P':
			# construct the name of the RBV file of the current residue
		        RBVFile=args.basisSet+"/Ac_"+res+"P_NHMe/rEV_norm.dat"
		else:
		        # construct the name of the RBV file of the current residue
		        RBVFile=args.basisSet+"/Ac_"+res+"_NHMe/rEV_norm.dat"
		print res+": "+RBVFile        
		# read RBV file
		with open(RBVFile, 'r') as f:
        		# add RBVs to the list of RBVs 
       			RBV.append(np.loadtxt(f))
	else:
	        # a zero is added to the RBV list if the residue is not included
	        print res+": excluded"       
	        RBV.append(np.array([0]))
#convert to numpy array
RBV=np.array(RBV)           
print ("---------------------------------------------------")
     

#------------------------------------------------
# Read microstate trajectories  
#------------------------------------------------
print "Reading microstate trajectories"
# List of microstate trajectories
# has to be a python list because the items have different lengths
# i.e. microstate traj for included residues, zero for excluded residues
microstateTraj=[]
# list of traj length
trajLengthList=[]
# loop over residues
# the microstate traj of each residue are added to the list
for i, res in enumerate(sequence):
    if resIncluded[i]==1:
        # construct the name of the RBV file of the current residue
        microstateTrajFile=args.inDir+"/"+str(i+1)+res+".npy"
        print res+": "+microstateTrajFile        
        # read microstate traj file
        with open(microstateTrajFile, 'r') as f:
           # add RBVs to the list of RBVs 
           thisTraj = np.load(f)
           # add current traj length to list of traj lengths
           trajLengthList.append(thisTraj.size)
           # check whether the data type of the trajectory is integer:
           if np.issubclass_(thisTraj.dtype.type, np.integer):
               microstateTraj.append(thisTraj)
           else:
               print "WARNING: converting data type to integer"
               microstateTraj.append(thisTraj.astype(int))             
    else:
       # a zero is added to the RBV list if the residue is not included
       print res+": excluded"        
       microstateTraj.append([0])
#convert to numpy array
microstateTraj=np.array(microstateTraj) 

# convert to numpy array
trajLengthList=np.array(trajLengthList)
# get the shortest trajectory length
trajLength=np.min(trajLengthList)
# check if the list contains trajectories with different lengths
if np.any(trajLengthList!=trajLength):
    print "WARNING: Trajectories have different lengths. Using the shortest trajectory length."
print ("---------------------------------------------------")

#------------------------------------------------
# Convert microstate traj to RBV traj  
#------------------------------------------------
print "Constructing RBV trajectories"

PBVTraj=np.ones(trajLength)    
    
# loop over positions in the current pbv
for r, thisRes in enumerate(sequence):
	# if the residue is included
        if resIncluded[r]==1:

		# define highest excitation level
		if thisRes=='P':
			maxEx=2
		else:
			maxEx=3 

		# loop over all excitations
		for ex in range(maxEx):
  	         	try: 
        	       		thisRBV=RBV[r][:,ex]
        	   	except IndexError: 
        	       		print "ERROR: Higher excitation required than available with this basis set. Aborting."
        	       		sys.exit(1)

          		# project the microstate traj of residue r onto thisRBV
	   		try:
               			thisRBVTraj=thisRBV[microstateTraj[r]]
	   		except IndexError:
	       			print 'ERROR: Number of microstates and basis vector dimension do not match. Aborting.'
	       			sys.exit(1)

  			# Construct the name of the output file
    			RBVFile = args.outDir+'/'+str(r+1)+'_RBV'+str(ex+1)+'.npy' 
    			# write to file
    			np.save(RBVFile, thisRBVTraj)
    
print ("---------------------------------------------------")        
print "END: "+str(datetime.now())
print ("---------------------------------------------------")
print ""
