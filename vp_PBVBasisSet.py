# -*- coding: utf-8 -*-
"""
Created on Sun May 17 20:25:36 2015

@author: bettina
"""

#------------------------------------------------
#  Imports
#------------------------------------------------
import argparse
import numpy as np
#import os
import sys
import itertools
from datetime import datetime

#--------------------------------
__author__ = 'bettina'

# output to screen  
print ("---------------------------------------------------")
print ("vp_generatePBV.py")
print ("Generates indices for peptide-based vectors (PBVs)")
print ("B. Keller (bettina.keller@fu-berlin.de), May. 2015")
print ("---------------------------------------------------") 
print "START: "+str(datetime.now())
print ("---------------------------------------------------")

#------------------------------------------------
# Command line arguments  
#------------------------------------------------
parser = argparse.ArgumentParser(description='Generates indices for peptide-based vectors (PBVs)')
parser.add_argument('-outDir','--outDir', help='Output directory',required=True)
parser.add_argument('-sequence','--sequence', help='File with peptide/protein sequence and included residues',required=True)
parser.add_argument('-proline', action='store_true', default=False, help='Only use the first excitation level for proline. Default: False')
parser.add_argument('-single', action='store_true', default=False, help='Singly excited states. Default: False')
parser.add_argument('-double', action='store_true', default=False, help='Doubly excited states. Default: False')
parser.add_argument('-triple', action='store_true', default=False, help='Triply excited states. Default: False')
parser.add_argument('-quadruple', action='store_true', default=False, help='Fourfold excited states. Default: False')
parser.add_argument('-quintuple', action='store_true', default=False, help='Fivefold excited states. Default: False')
args = parser.parse_args()


print ("Using the following command line arguments")
print ("-outDir (Output directory):")
print ("%s" % args.outDir)
print ("-sequence (File with peptide/protein sequence and included residues):")
print ("%s" % args.sequence)
print ("-proline (Only use the first excitation level for proline. Default: False):")
print ("%s" % args.proline)
print ("-single (Singly excited states. Default: False):")
print ("%s" % args.single)
print ("-double (Doubly excited states. Default: False):")
print ("%s" % args.double)
print ("-triple (Triply excited states. Default: False):")
print ("%s" % args.triple)
print ("-quadruple (Fourfold excited states. Default: False):")
print ("%s" % args.quadruple)
print ("-quintuple (Fivefold excited states. Default: False):")
print ("%s" % args.quintuple)

print ("---------------------------------------------------")

#------------------------------------------------
# Functions   
#------------------------------------------------
#---------------
# Yields a list of tuples of length nExcitations 
# containing the various combinations of excitation levels.
# The list returned as a numpy-array
def getExcitationTuples( excitationLevels, nExcitations ):
   iterables = [ excitationLevels ]
   excitationTuples=[]
   for t in itertools.product(*iterables, repeat=nExcitations):
       excitationTuples.extend([t]);    
#   print x    
   return np.array(excitationTuples);

#---------------
# Yields a list of tuples of length nExcitations 
# containing the indices of the excited residues.
# The list returned as a numpy-array
def getExcitatedRes( nResIncluded, nExcitations ):
   resList = range(nResIncluded)
   resExcited=[]
   for t in itertools.combinations(resList, nExcitations):
       resExcited.extend([t]);    
#   print x    
   return np.array(resExcited);

#---------------
# Generates a PBV for each combination of
# tuple of excited residues and
# combination of excitation levels
def getPBV( excitedRes, excitationTuples, nResIncluded) :
    PBV=[]
    #loop over all tuples of excited residues
    for thisExcitedRes in excitedRes:
        # loop over all combinations of excitation levels
        for thisExcitationTuple in excitationTuples:
            # generate a ground state PBV 
            # as starting point for thisPBV
            thisPBV=np.empty(nResIncluded)
            thisPBV.fill(1)
            # loop over residues in thisExcitedRes
            # fill in with the excitation of thisExcitationTuple
            for index, r in enumerate(thisExcitedRes):
                thisPBV[r] = thisExcitationTuple[index]         
            # add thisPBV to the list of PBVs   
            PBV.append(thisPBV)        
            #print(thisPBV)
    return np.array(PBV);        

#---------------
# Converts PBVs for only included residues to PBVs for the full peptide
# by adding columns of zeros at the positions of the excluded residues
def addZerosForExcludedRes( PBV, resIncluded) :
    # get the position/index of the included residues
    resIncludedPos=np.where(resIncluded==1)[0]
    # total number of residues
    nRes=resIncluded.size

    # only one PBV in the list of PBVs
    if PBV.ndim==1:
        # check whether the length of the PBVs match the number of included res.
        if resIncludedPos.size!=PBV.size:
             print "ERROR: PBV does not match number of included residues."
             sys.exit(1)
        # number of PBV
        nPBV=1
        # initialize new (full) PBV array
        PBV_full=np.zeros(nRes)
        # loop over all res
        for counter, pos in enumerate(resIncludedPos): 
            PBV_full[pos]=PBV[counter]                    
    # more than one PBV    
    elif PBV.ndim==2:
        # check whether the length of the PBVs match the number of included res.
        if resIncludedPos.size!=PBV.shape[1]:
             print "ERROR: PBV does not match number of included residues."
             sys.exit(1)
        # number of PBV
        nPBV=PBV.shape[0]
        # initialize new (full) PBV array
        PBV_full=np.zeros(nPBV*nRes).reshape(nPBV,nRes)
        # loop over all res
        for counter, pos in enumerate(resIncludedPos): 
            PBV_full[:,pos]=PBV[:,counter]
            
    # return the full PBV
    return PBV_full;    

#---------------
# Deletes PBVs from the list of PBVs in which excitations higher than 3
# are assigned to a proline residue
def deleteExcitationsForPro( PBV, sequence) :
    # extract the actual sequence from the array sequence
    sequence=np.array(zip(*sequence)[0])
    # get the positions of the prolines
    posPro=np.where(sequence=='P')[0]
    print "Proline: deleting excitations > 2 for residue(s) "+str(posPro+1)
    # loop over prolines
    for pos in posPro:
        # get row numbers in which excitations are higher than 2
        rows=np.where(PBV[:,pos]>2)
        PBV=np.delete(PBV,rows,0)
    return PBV;

#------------------------------------------------
# excitated states per residue   
#------------------------------------------------
# hardcoded to [2,3] 
# a switch via a command line argument to add higher excitation levels
# should be implemented here
excitationLevels=[2,3]

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


#------------------------------------------------
# Generate excited states  
#------------------------------------------------
# ground state
PBV00=np.empty(nResIncluded)
PBV00.fill(1)
# add a column of zeros for the excluded residues
PBV00=addZerosForExcludedRes( PBV00, resIncluded) 

# singly excited states
if args.single==True and nResIncluded >=1:
    print("Generating singly excited states.")
    # generate all tuples of excitatios of length 1
    excitationTuples01 =  getExcitationTuples( excitationLevels, 1 ) 
    # generate all combinations of residues indeces of length 1
    #(w/o replacment) 
    excitedRes01 = getExcitatedRes( nResIncluded, 1)
    # combine all excitation tuples with all residue combinatons
    # to get the corresponding peptide-based vectors
    PBV01 = getPBV(excitedRes01, excitationTuples01, nResIncluded)
    # add a column of zeros for the excluded residues
    PBV01=addZerosForExcludedRes( PBV01, resIncluded) 
    # delete higher excited states for prolines
    if args.proline==True:
        PBV01=deleteExcitationsForPro( PBV01, sequence) 

elif args.single==True and nResIncluded <1:    
    print "ERROR: singly excited states"
    print "ERROR: number of included residues smaller than number of excited residues"
    sys.exit(1)

# doubly excited states    
if args.double==True and nResIncluded >=2:
    print("Generating doubly excited states.")
    # generate all tuples of excitatios of length 2
    excitationTuples02 =  getExcitationTuples( excitationLevels, 2 ) 
    # generate all combinations of residues indeces of length 2
    #(w/o replacment) 
    excitedRes02 = getExcitatedRes( nResIncluded, 2)
    # combine all excitation tuples with all residue combinatons
    # to get the corresponding peptide-based vectors
    PBV02 = getPBV(excitedRes02, excitationTuples02, nResIncluded)    
    # add a column of zeros for the excluded residues
    PBV02=addZerosForExcludedRes( PBV02, resIncluded) 
    # delete higher excited states for prolines
    if args.proline==True:
        PBV02=deleteExcitationsForPro( PBV02, sequence) 

elif args.double==True and nResIncluded <2:    
    print "ERROR: doubly excited states"
    print "ERROR: number of included residues smaller than number of excited residues"
    sys.exit(1)

# triply excited states
if args.triple==True and nResIncluded >=3:
    print("Generating triply excited states.")
    # generate all tuples of excitatios of length 3
    excitationTuples03 =  getExcitationTuples( excitationLevels, 3 ) 
    # generate all combinations of residues indeces of length 3
    #(w/o replacment) 
    excitedRes03 = getExcitatedRes( nResIncluded, 3)    
    # combine all excitation tuples with all residue combinatons
    # to get the corresponding peptide-based vectors
    PBV03 = getPBV(excitedRes03, excitationTuples03, nResIncluded)    
    # add a column of zeros for the excluded residues
    PBV03=addZerosForExcludedRes( PBV03, resIncluded) 
    # delete higher excited states for prolines
    if args.proline==True:
        PBV03=deleteExcitationsForPro( PBV03, sequence) 
        

elif args.triple==True and nResIncluded <3:    
    print "ERROR: triply excited states"
    print "ERROR: number of included residues smaller than number of excited residues"
    sys.exit(1)

# fourfold excited states
if args.quadruple==True and nResIncluded >=4:
    print("Generating fourfold excited states.")
    # generate all tuples of excitatios of length 4
    excitationTuples04 =  getExcitationTuples( excitationLevels, 4 ) 
    # generate all combinations of residues indeces of length 4
    # (w/o replacment) 
    excitedRes04 = getExcitatedRes( nResIncluded, 4)        
    # combine all excitation tuples with all residue combinatons
    # to get the corresponding peptide-based vectors
    PBV04 = getPBV(excitedRes04, excitationTuples04, nResIncluded)    
    # add a column of zeros for the excluded residues
    PBV04=addZerosForExcludedRes( PBV04, resIncluded) 
    # delete higher excited states for prolines
    if args.proline==True:
        PBV04=deleteExcitationsForPro( PBV04, sequence) 

elif args.quadruple==True and nResIncluded <4:    
    print "ERROR: fourfold excited states"
    print "ERROR: number of included residues smaller than number of excited residues"
    sys.exit(1)

# fivefold excited states
if args.quintuple==True and nResIncluded >=5:
    print("Generating fivefold excited states.")
    # generate all tuples of excitatios of length 5
    excitationTuples05 =  getExcitationTuples( excitationLevels, 5 ) 
    # generate all combinations of residues indeces of length 5
    # (w/o replacment) 
    excitedRes05 = getExcitatedRes( nResIncluded, 5)  
    # combine all excitation tuples with all residue combinatons
    # to get the corresponding peptide-based vectors
    PBV05 = getPBV(excitedRes05, excitationTuples05, nResIncluded)          
     # add a column of zeros for the excluded residues
    PBV05=addZerosForExcludedRes( PBV05, resIncluded) 
    # delete higher excited states for prolines
    if args.proline==True:
        PBV05=deleteExcitationsForPro( PBV05, sequence) 

elif args.quintuple==True and nResIncluded <5:    
    print "ERROR: fivefold excited states"
    print "ERROR: number of included residues smaller than number of excited residues"
    sys.exit(1)


#------------------------------------------------
# Write PBVs to file
#------------------------------------------------
# Generate output file name
outFile=['x','x','x','x','x']
if args.single==True:
    outFile[0]='S'
if args.double==True:
    outFile[1]='D'
if args.triple==True:
    outFile[2]='T'
if args.quadruple==True:
    outFile[3]='Q'
if args.quintuple==True:
    outFile[4]='Q'
outFile=args.outDir+"/PBV_"+"".join(outFile)+".dat"    

with open(outFile, 'w') as f:
    # write the ground state to file    
    #np.savetxt(f,PBV00[None],delimiter="\t", fmt='%i', header="ground state")
    np.savetxt(f,PBV00[None],delimiter="\t", fmt='%i')
    # singly excited state
    if args.single==True:
        #np.savetxt(f,PBV01,delimiter="\t", fmt='%i', header="singly excited states")
	np.savetxt(f,PBV01,delimiter="\t", fmt='%i')
    # doubly excited state
    if args.double==True:
        #np.savetxt(f,PBV02,delimiter="\t", fmt='%i', header="doubly excited states")
	np.savetxt(f,PBV02,delimiter="\t", fmt='%i')
    # triply excited state
    if args.triple==True:
        #np.savetxt(f,PBV03,delimiter="\t", fmt='%i', header="triply excited states")
	np.savetxt(f,PBV03,delimiter="\t", fmt='%i')
    # fourfold excited state
    if args.quadruple==True:
        #np.savetxt(f,PBV04,delimiter="\t", fmt='%i', header="quadruply excited states")
	np.savetxt(f,PBV04,delimiter="\t", fmt='%i')
    # fivefold excited state
    if args.quintuple==True:
        #np.savetxt(f,PBV05,delimiter="\t", fmt='%i', header="fivefold excited states")
	np.savetxt(f,PBV05,delimiter="\t", fmt='%i')
print ("---------------------------------------------------")        
print "END: "+str(datetime.now())
print ("---------------------------------------------------")
print ""
