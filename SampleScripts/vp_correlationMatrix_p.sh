#!/bin/sh

# Variables
SYSTEM=hIAPP
PBVSET=SD
PBVMOD=PBV
JOBS=10

# Generate directory names
WORKDIR=/home/louism/Documents/${SYSTEM}
PYTHONDIR=/home/louism/Documents/VariationalApproach/variational

# pbvfile
PBVFILE=PBV_SDxxx.dat

# overlap matrix
python  $PYTHONDIR/vp_correlationMatrix_p.py \
		-inDir	${WORKDIR}/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep1 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep2 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep3 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep4 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep5 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep6 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep7 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep8  $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep9 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep10 \
		-outDir $WORKDIR/correlationMatrix/${PBVMOD}_${PBVSET} \
		-PBVFile $WORKDIR/pbvsets/${PBVFILE} \
		-overlapMatrix \
		-nJobs ${JOBS}

# correlation matrix
for TAU in 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 20000 30000 40000 50000 60000 70000 80000 90000 100000; do
        python $PYTHONDIR/vp_correlationMatrix_p.py \
		-inDir	${WORKDIR}/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep1 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep2 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep3 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep4 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep5 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep6 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep7 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep8  $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep9 $WORKDIR/nobackup/PBVTrace/${PBVMOD}_${PBVSET}/rep10 \
		-outDir $WORKDIR/correlationMatrix/${PBVMOD}_${PBVSET} \
           	-PBVFile $WORKDIR/pbvsets/${PBVFILE} \
		-reversible \
		-lagTime $TAU \
		-nJobs ${JOBS}
done

