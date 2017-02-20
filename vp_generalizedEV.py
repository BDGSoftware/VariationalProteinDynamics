# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 11:36:55 2015

@author: bettina
"""



#------------------------------------------------
#  Imports
#------------------------------------------------
import argparse
import numpy as np
from scipy import linalg
#import os
import sys
from datetime import datetime

#--------------------------------
__author__ = 'bettina'

# output to screen  
print ("---------------------------------------------------")
print ("vp_correlationMatrix.py")
print ("Solves the generalized eigenvalue problem for a correlation and an overlap matrix")
print ("B. Keller (bettina.keller@fu-berlin.de), August 2015")
print ("---------------------------------------------------") 
print "START: "+str(datetime.now())
print ("---------------------------------------------------")

#------------------------------------------------
# Command line arguments  
#------------------------------------------------
parser = argparse.ArgumentParser(description='Solves the generalized eigenvalue problem for a correlation and an overlap matrix')
parser.add_argument('-S','--S', help='File which contains the overlap matrix (npy)',required=True)
parser.add_argument('-C','--C', help='File which contains the correlation matrix (npy)',required=True)
parser.add_argument('-eigenValueFile','--eigenValueFile', help='File to which the eigenvalues are written',required=True)
parser.add_argument('-eigenVectorFile','--eigenVectorFile', help='File to which the eigenvectors are written',required=True)
parser.add_argument('-ascii', action='store_true', default=False, help='Write output to ascii-files. Default: False')


args = parser.parse_args()

print ("Using the following command line arguments")
print ("-S (File which contains the overlap matrix (npy)):")
print ("%s" % args.S) 
print ("-C (File which contains the correlation matrix (npy)):")
print ("%s" % args.C)
print ("-eigenValueFile (File to which the eigenvalues are written):")
print ("%s" % args.eigenValueFile)
print ("-eigenValueVector (File to which the eigenvectors are written):")
print ("%s" % args.eigenVectorFile)
print ("-ascii (Write output to ascii-files. Default: False):")
print ("%s" % args.ascii)


print ("---------------------------------------------------")


#------------------------------------------------
# Read overlap matrix and correlation matrix  
#------------------------------------------------
with open(args.S, 'r') as f:
        S = np.load(f)
        
with open(args.C, 'r') as f:
        C = np.load(f)
        
        
#------------------------------------------------
# Sanity checks for the matrices  
#------------------------------------------------        
if S.shape!=C.shape: 
    print "ERROR: Overlap matrices and correlation matrix do not have the same dimension."
    sys.exit(1)

# check if we have symmetric matrices
symmetricS=False
if np.array_equal(S, S.transpose()):
    symmetricS=True

symmetricC=False
if np.array_equal(C, C.transpose()):
    symmetricC=True

# Check whether the overlap matrix has full rank    
if S.shape[1]!=np.linalg.matrix_rank(S):
    rankdef=True
    print "WARNING: Overlap matrix does not have full rank."
else:
    rankdef=False

# check overlap matrix for positive definiteness
try:
    linalg.cholesky(S)
    pd=True
except linalg.LinAlgError:
    pd=False
    print 'WARNING: Overlap matrix is not positive definite.'
    
#------------------------------------------------
# Solve generalized eigenvalue problem
#------------------------------------------------     
if symmetricC==True and symmetricS==True:
    print 'Correlation matrix and overlap matrix symmetric.'
    if rankdef or not pd:
        from pyemma.util.linalg import eig_corr
	print 'Using pyemma.util.linalg.eig_corr()'
	eigenValues,eigenVectors = eig_corr(S,C)
    else:
        print "Using scipy.linalg.eig()"
    	eigenValues,eigenVectors = linalg.eigh(C,S)
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]        
else:
    print "Correlation matrix and/or overlap matrix not symmetric."
    if rankdef or not pd:
	print 'Using pyemma.util.linalg.eig_corr()'
	eigenValues,eigenVectors = eig_corr(S,C)
    else:
        print "Using scipy.linalg.eig()"
    	eigenValues,eigenVectors = linalg.eigh(C,S)
    idx = eigenValues.argsort()[::-1]   
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx] 
    
    
#------------------------------------------------
# Write matrix to file
#------------------------------------------------
print "Writing eigenvalues to: "
print args.eigenValueFile
if(args.ascii==False): 
    np.save(args.eigenValueFile, eigenValues)
else: 
    np.savetxt(args.eigenValueFile, eigenValues)

print "Writing eigenvectors to: "
print args.eigenVectorFile
if(args.ascii==False): 
    np.save(args.eigenVectorFile, eigenValues)
else: 
    np.savetxt(args.eigenVectorFile, eigenVectors)

print ("---------------------------------------------------")
print "END: "+str(datetime.now())
print ("---------------------------------------------------")
print ""

    
