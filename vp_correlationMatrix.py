# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 11:30:29 2015

@author: bettina
"""


#------------------------------------------------
#  Imports
#------------------------------------------------
import argparse
import numpy as np
#import os
import sys
from datetime import datetime

#--------------------------------
__author__ = 'bettina'

# output to screen  
print ("---------------------------------------------------")
print ("vp_correlationMatrix.py")
print ("Calculates the time-lagged correlation matrix or the overlap matrix of the a set of PBV trajectories")
print ("B. Keller (bettina.keller@fu-berlin.de), July. 2015")
print ("---------------------------------------------------") 
print "START: "+str(datetime.now())
print ("---------------------------------------------------")

#------------------------------------------------
# Command line arguments  
#------------------------------------------------
parser = argparse.ArgumentParser(description='Calculates the time-lagged correlation matrix of the a set of PBV trajectories')
parser.add_argument('-inDir','--inDir', nargs='*', help='List of directories with PBV-trajectories (npy)',required=True)
parser.add_argument('-outDir','--outDir', help='Directory to which the overlap/correlation matrix is written',required=True)
parser.add_argument('-PBVFile','--PBVFile', help='File with PBV indices',required=True)
#
parser.add_argument('-overlapMatrix', action='store_true', default=False, help='Calculate the overlap matrix instead of a correlation matrix. Default: False')
parser.add_argument('-reversible', action='store_true', default=False, help='Enforce reversibility by symmetrizing the matrix. Can be ommited / is ignored if overlapMatrix==True.Default: False')
parser.add_argument('-lagTime','--lagTime', type=int, help='Lag time in time steps of the PBV trajectory. Can be ommitted / is ignored if overlapMatrix==True',required=False)
parser.add_argument('-PBVBlock','--PBVBlock',type=int, help='Use blocks of PBV trajectories to construct the overlap/correlation matrix. PBVBlock > 1 can speed up the program. Default: PBVBlock=1',required=False)


args = parser.parse_args()


print ("Using the following command line arguments")
print ("-inDir (List of directories with PBV-trajectories (npy)):")
print ("%s" % args.inDir) 
print ("-outDir (Directory to which the overlap/correlation matrix is written):")
print ("%s" % args.outDir)
print ("-PBVFile (File with PBV indices):")
print ("%s" % args.PBVFile)
print ("-overlapMatrix (Calculate the overlap matrix instead of a correlation matrix. Default: False):")
print ("%s" % args.overlapMatrix)
print ("-reversible (Enforce reversibility by symmetrizing the matrix.  Can be ommited / is ignored if overlapMatrix==True. Default: False):")
print ("%s" % args.reversible)
print ("-lagTime (Lag time in time steps of the PBV trajectory. Can be ommitted / is ignored if overlapMatrix==True):")
print ("%s" % args.lagTime)
print ("-PBVBlock (Use blocks of PBV trajectories to construct the overlap/correlation matrix. PBVBlock > 1 can speed up the program. Default: PBVBlock=1):")
print ("%s" % args.PBVBlock)

print ("---------------------------------------------------")

#------------------------------------------------
# Sanity checks 
#------------------------------------------------
if args.lagTime!=None and args.lagTime<0: 
    print "ERROR: lag time negative"
    sys.exit(1)

if args.PBVBlock!=None: 
    print "WARNING: PBV block not yet implemented. Proceeding with PBVBlock=1"

if args.PBVBlock!=None and args.PBVBlock<0: 
    print "ERROR: PBV block negative"
    sys.exit(1)

if args.overlapMatrix==True and args.lagTime==None:
    print "Estimating an overlap matrix." 
elif args.overlapMatrix==True and args.lagTime!=None:    
    print "Estimating an overlap matrix. Ignoring the lag time information."
elif args.overlapMatrix==False and args.lagTime==None:
    print "ERROR: No lag time specified for the estimation of the correlation matrix"
    sys.exit(1)
elif args.overlapMatrix==False and args.lagTime!=None:
    print ("Estimating a correlation matrix with lag time: %s" % args.lagTime)
          
print ("---------------------------------------------------")

#------------------------------------------------
# Read PBV index file  
#------------------------------------------------
PBVIndices=[]
# open file
with open(args.PBVFile, 'r') as f: 
    #loop over lines
    for line in f:
        # check wether the line is a comment        
        if not line.startswith("#"):
            PBVIndices.append(line.strip().split("\t")) 


#------------------------------------------------
# Parameters
#------------------------------------------------
# set lag time
if args.overlapMatrix==True: 
    tau=0
else:   
    tau=args.lagTime 

# number of PBV trajectories
nPBV=len(PBVIndices)


#------------------------------------------------
# Estimate
#------------------------------------------------
if(args.overlapMatrix==False): 
    print "Reading PBV trajectories and estimating correlation matrix."
else:
    print "Reading PBV trajectories and estimating overlap matrix."    
#list of correlation matrices
CList=[]
#list of trajectory lengths
lTrajList=[]

# loop over PBV directories and estimate the correlation matrix
# for each trajectoriy
for d, PBVDir in enumerate(args.inDir):    
    print(PBVDir)    
    
    # initialize overlap/correlation matrix
    thisC=np.zeros((nPBV,nPBV))
    
    # load sample trajectory length
    with open(PBVDir+"/PBV"+"".join(PBVIndices[0])+".npy", 'r') as f:
        pbv = np.load(f)
    # get trajectory length    
    lTraj=pbv.size     
    
    # if trajectory is too short
    if(lTraj<=tau):
        # issue a warning and skip this set of PBV trajectories
        print "WARNING: PBV trajectories in " + PBVDir + " shorter than lag time. Skipping these trajectories."    

    # if trajectory is long enough continue    
    else:    
        # add this trajectory length to the list of trajectory lengths
        lTrajList.append(lTraj)
        # delete sample trajectory
        del pbv
                
        # loop over pairs of PBV trajectories    
        for i in range(len(PBVIndices)):
            # read PBV traj file i
            with open(PBVDir+"/PBV"+"".join(PBVIndices[i])+".npy", 'r') as fi:
                pbvi = np.load(fi)
    
            for j in range(i, len(PBVIndices)):
                # read PBV traj file i
                with open(PBVDir+"/PBV"+"".join(PBVIndices[j])+".npy", 'r') as fj:
                    # add RBVs to the list of RBVs 
                    pbvj = np.load(fj)
                
                # calculate correlation
                # PBV traj i (0 to T-tau), PBV traj j (tau+1 to T)
                thisC[i,j]=np.dot(pbvi[0:lTraj-tau],pbvj[tau:lTraj])/lTraj
                
                # transposed element for correlation matrices            
                if args.overlapMatrix==False and i!=j:
                    thisC[j,i]=np.dot(pbvj[0:lTraj-tau],pbvi[tau:lTraj])/lTraj
                # transposed element for overlap matrices                            
                elif args.overlapMatrix==True and i!=j:     
                    thisC[j,i]=thisC[i,j]
                

    # add this correlation matrix to the list of correlation matrices
    CList.append(thisC)

# combine the correlation matrices
if(sum(lTrajList)>0):
    C=np.zeros([nPBV,nPBV])
    for thisC, thisTrajLength in zip(CList, lTrajList):
        C+=thisC*thisTrajLength
    C/=sum(lTrajList)
else: 
   print ("---------------------------------------------------")     
   print "ERROR: None of the PBV trajectories where longer than the lag time. No correlation matrix estimated." 
   sys.exit(1)
  
    
#------------------------------------------------
# Reversibility
#------------------------------------------------
# enforce reversibility for correlation matrices
# (overlap matrices are symmetric by construction)   
if args.overlapMatrix==False and args.reversible==True: 
    print ("---------------------------------------------------")    
    print "Enforcing reversibility"
    C=(C+C.transpose())/2


#------------------------------------------------
# Write matrix to file
#------------------------------------------------
# Generate output file name
if args.overlapMatrix==True:
    outFile=args.outDir+"/S.npy" 
else:
    if tau<10:
        outFile=args.outDir+"/C_tau00000"+str(tau)+".npy" 
    elif 10<=tau<100:
        outFile=args.outDir+"/C_tau0000"+str(tau)+".npy" 
    elif 10<=tau<1000:
        outFile=args.outDir+"/C_tau000"+str(tau)+".npy" 
    elif 10<=tau<10000:
        outFile=args.outDir+"/C_tau00"+str(tau)+".npy" 
    elif 10<=tau<100000:
        outFile=args.outDir+"/C_tau0"+str(tau)+".npy" 
    else:
        outFile=args.outDir+"/C_tau"+str(tau)+".npy" 

print ("---------------------------------------------------")    
if(args.overlapMatrix==False): 
    print "Writing correlation matrix to: "
    print outFile
else:
    print "Writing overlap matrix to: "
    print outFile
np.save(outFile, C)


print ("---------------------------------------------------")
print "END: "+str(datetime.now())
print ("---------------------------------------------------")
print ""

