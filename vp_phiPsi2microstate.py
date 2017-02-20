# -*- coding: utf-8 -*-
"""
Created on Sun May 17 11:58:37 2015

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
print ("vp_phiPsi2microstate.py")
print ("Projects phi-psi-trajectories onto a regular grid of microstates")
print ("B. Keller (bettina.keller@fu-berlin.de), May. 2015")
print ("---------------------------------------------------") 
print "START: "+str(datetime.now())
print ("---------------------------------------------------")

#------------------------------------------------
# Command line arguments  
#------------------------------------------------
parser = argparse.ArgumentParser(description='Projects phi-psi-trajectories onto a regular grid of microstates')
parser.add_argument('-inDir','--inDir', help='Directory with phi-psi-trajectories (ascii)',required=True)
parser.add_argument('-outDir','--outDir', help='Output directory',required=True)
parser.add_argument('-sequence','--sequence', help='File with peptide/protein sequence and included residues',required=True)
parser.add_argument('-shift', action='store_true', default=False, help='Shift the torsion angles from [-180,180[ to [0,360[ (Default:False)')
parser.add_argument('-tolerance', '--tolerance', default=0.0, type=float, help='Tolerance for deviations from [0,360[, default: 0.0',required=False)
parser.add_argument('-nBins','--nBins', help='Number of bins per torsion angle',required=True)
args = parser.parse_args()


print ("Using the following command line arguments")
print ("-inDir (Directory with phi-psi-trajectories (ascii)):")
print ("%s" % args.inDir) 
print ("-outDir (Output directory):")
print ("%s" % args.outDir)
print ("-sequence (File with peptide/protein sequence and included residues):")
print ("%s" % args.sequence)
print ("-shift (Shift the torsion angles from [-180,180[ to [0,360[ (Default: False)):")
print ("%s" % args.shift)
print ("-tolerance (Tolerance for deviations from [0,360[, default: 0.0:")
print ("%s" % args.tolerance)
print ("-nBins (Number of bins per torsion angle):")
print ("%s" % args.nBins)
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

#------------------------------------------------
# Create list of input and output file names  
#------------------------------------------------
inFileList=[]
outFileList=[]
for i, thisRes in enumerate(sequence): 
    if thisRes[1]=='1':
        inFileList.append(args.inDir+"/"+str(i+1)+thisRes[0]+".dat")
        outFileList.append(args.outDir+"/"+str(i+1)+thisRes[0]+".npy")
        
#print inFileList
#print outFileList


#------------------------------------------------
# Convert trajectories  
#------------------------------------------------
for i, thisFile in enumerate(inFileList):
    print thisFile
    
    # read phi-psi trajectory
    phiPsi = np.loadtxt(thisFile)
    
    # shift to angle values to [0,360[
    if args.shift==True:
        phiPsi+=180        
    
    # check whether tolerance is violated
    if np.any(phiPsi>360+args.tolerance) or np.any(phiPsi<=0-args.tolerance):
            print "ERROR: phi-psi values out of bound"
            sys.exit(1)

    # shift values which are within [0-tolerance, 360+tolerance[ to [0,360[    
    phiPsi=np.where(phiPsi>=360,359.9999,phiPsi)
    phiPsi=np.where(phiPsi<0,0,phiPsi)
    
    #set bin width
    binWidth=360/float(args.nBins)
    
    # convert to microstates
    microstates=np.floor(phiPsi[:,0]/binWidth)*float(args.nBins) + np.floor(phiPsi[:,1]/binWidth)
    
    # write to file
    np.save(outFileList[i], microstates.astype(int))


print ("---------------------------------------------------")        
print "END: "+str(datetime.now())
print ("---------------------------------------------------")
print ""
