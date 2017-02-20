# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 10:44:58 2014

@author: bettina
"""

#------------------------------------------------
#  Imports
#------------------------------------------------
import argparse
import numpy as np
import os
from datetime import datetime


#--------------------------------
__author__ = 'bettina'

# output to screen  
print ("---------------------------------------------------")
print ("vp_npy2dat.py")
print ("Converts a list of npy-files to ascii-files")
print ("B. Keller (bettina.keller@fu-berlin.de), August 2015")
print ("---------------------------------------------------") 
print "START: "+str(datetime.now())
print ("---------------------------------------------------")

#------------------------------------------------
# Command line arguments  
#------------------------------------------------
parser = argparse.ArgumentParser(description='Converts file from npy to ascii. B Keller, August 2015')
parser.add_argument('-inDir','--inDir', help='Directory with ascii files',required=True)
parser.add_argument('-outDir','--outDir', help='Output directory',required=True)
parser.add_argument('-fileList','--fileList', help='List with ascii files',required=True)
args = parser.parse_args()


print ("Using the following command line arguments")
print ("-inDir (Directory with ascii files):")
print ("%s" % args.inDir) 
#print
print ("-outDir (Output directory):")
print ("%s" % args.outDir)
#print
print ("-fileList (List with ascii files):")
print ("%s" % args.fileList)
print ("---------------------------------------------------")


#------------------------------------------------
# Read files and convert to npy-format  
#------------------------------------------------
# open file
with open(args.fileList, 'r') as f: 
    #loop over lines
    for line in f:
       # strip new-line sign
       # construct full input file name            
       inFile = args.inDir+"/"+line.strip()
       # read npy file
       x = np.load(inFile)
       # construct full output file name
       outFile = args.outDir +"/"+os.path.splitext(line.strip())[0]+".dat"
       # save npy file
       np.savetxt(outFile, x)
       print outFile

print ("---------------------------------------------------")        
print "END: "+str(datetime.now())
print ("---------------------------------------------------")
print ""