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
print ("vp_dat2npy.py")
print ("Converts a list of ascii-files to npy-files")
print ("B. Keller (bettina.keller@fu-berlin.de), May. 2015")
print ("---------------------------------------------------") 
print "START: "+str(datetime.now())
print ("---------------------------------------------------")

#------------------------------------------------
# Command line arguments  
#------------------------------------------------
parser = argparse.ArgumentParser(description='Converts file from ascii to npy. B Keller, May 2015')
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
        # check wether the line is a comment        
        if not line.startswith("#"):
            # strip new-line sign
            # construct full input file name            
            inFile = args.inDir+"/"+line.strip()
            # read ascii file
            x = np.loadtxt(inFile)
            # construct full output file name
            outFile = args.outDir +"/"+os.path.splitext(line.strip())[0]+".npy"
            # save npy file
            np.save(outFile, x)
            print outFile

print ("---------------------------------------------------")        
print "END: "+str(datetime.now())
print ("---------------------------------------------------")
print ""
