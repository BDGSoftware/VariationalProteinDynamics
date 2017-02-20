
#------------------------------------------------
#  Imports
#------------------------------------------------
import argparse
import numpy as np
#import os
import sys
from datetime import datetime

# output to screen  
print ("---------------------------------------------------")
print ("vp_phiPsi23state.py")
print ("Projects phi-psi-trajectories onto a user-defined grid of microstates")
print ("---------------------------------------------------") 
print "START: "+str(datetime.now())
print ("---------------------------------------------------")

#------------------------------------------------
# Command line arguments  
#------------------------------------------------
parser = argparse.ArgumentParser(description='Projects phi-psi-trajectories onto a user-defined grid of microstates')
parser.add_argument('-inDir','--inDir', help='Directory with phi-psi-trajectories (ascii)',required=True)
parser.add_argument('-outDir','--outDir', help='Output directory',required=True)
parser.add_argument('-sequence','--sequence', help='File with peptide sequence, included residues and coordinates of boundary bins (format: phi_left, phi_right, psi_bottom, psi_top)',required=True)
parser.add_argument('-shift', action='store_true', default=False, help='Shift the torsion angles from [-180,180[ to [0,360[ (Default:False)')
parser.add_argument('-tolerance', '--tolerance', default=0.0, type=float, help='Tolerance for deviations from [0,360[, default: 0.0',required=False)
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
print ("---------------------------------------------------")

#------------------------------------------------
# Binning algorithms
#------------------------------------------------
def mstate(trj,coords):    
        c_phi=2*np.multiply((trj[:,0]>=coords[0]),(trj[:,0]<coords[1]))
        c_psi=1*(trj[:,1]<=coords[2])+1*(trj[:,1]>coords[3])
	c_com=c_phi+c_psi
        c_com=np.where(c_com==3,2,c_com)
        return c_com

def mstate_p(trj,coords):
	c_com=1*(trj[:,1]<=coords[0])+1*(trj[:,1]>coords[1])
	return c_com

def mstate_g(trj,coords):
	if coords.shape[0]==5:
		coords=np.hstack((coords,np.array([coords[4],coords[5]])))
	elif not coords.shape[0]==7:
		print "ERROR: Wrong format of glycine binning coordinates"
            	sys.exit(1)

	phi_l=-2*np.multiply((trj[:,0]>=coords[0]),(trj[:,0]<coords[1]))
	phi_r=-1*np.multiply((trj[:,0]>=coords[1]),(trj[:,0]<coords[2]))
	psi_mo=4*np.multiply((trj[:,1]>=coords[3]),(trj[:,1]<coords[4]))*(phi_r+1)
	psi_mi=-4*np.multiply((trj[:,1]>=coords[5]),(trj[:,1]<coords[6]))*phi_r
	psi=psi_mo+psi_mi
	c_com=phi_l+phi_r+psi
	c_com=np.where(c_com==0,5,c_com)
	c_com=np.where(c_com==-1,1,c_com)
	c_com=np.where(c_com==-2,0,c_com)
	return c_com
#------------------------------------------------
# Read sequence file  
#------------------------------------------------
sequence=[]
# open file
with open(args.sequence, 'r') as f: 
    #loop over lines
    for line in f:
        # check whether the line is a comment        
        if not line.startswith("#"):
            sequence.append(line.strip().split("\t"))

#------------------------------------------------
# Create list of input and output file names  
#------------------------------------------------
inFileList=[]
outFileList=[]
for i, thisRes in enumerate(sequence): 
    if thisRes[1]=='1':
        inFileList.append([args.inDir+"/"+str(i+1)+thisRes[0]+".dat", [int(thisRes[k]) for k in range(2,len(thisRes))]])
        outFileList.append(args.outDir+"/"+str(i+1)+thisRes[0]+".npy")
    else:
        sequence=np.delete(sequence,i)
        
#------------------------------------------------
# Convert trajectories  
#------------------------------------------------
for i, thisFile in enumerate(inFileList):
    print thisFile
    
    # read phi-psi trajectory and coordinates
    phiPsi = np.loadtxt(thisFile[0])
    coords = np.array(thisFile[1])

    # shift to angle values to [0,360[
    if args.shift==True:
        phiPsi+=180        
        if np.any(coords<0):
            print 'ERROR: coordinate values out of bound'
            sys.exit(1)
    elif np.any(coords>180):
            print 'ERROR: coordinate values out of bound'
            sys.exit(1)

    # check whether tolerance is violated
    if np.any(phiPsi>360+args.tolerance) or np.any(phiPsi<=0-args.tolerance):
            print "ERROR: phi-psi values out of bound"
            sys.exit(1)

    # shift values which are within [0-tolerance, 360+tolerance[ to [0,360[    
    phiPsi=np.where(phiPsi>=360,359.9999,phiPsi)
    phiPsi=np.where(phiPsi<0,0,phiPsi)
    
    # convert to microstates
    if sequence[i][0]=='P':
	microstates=mstate_p(phiPsi,coords)
 
    elif sequence[i][0]=='G':
    	microstates=mstate_g(phiPsi,coords)
    
    else:
	microstates=mstate(phiPsi,coords)
    
    # write to file
    np.save(outFileList[i], microstates)


print ("---------------------------------------------------")        
print "END: "+str(datetime.now())
print ("---------------------------------------------------")
print ""
