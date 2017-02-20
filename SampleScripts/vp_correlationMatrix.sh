#!/bin/sh

WORKDIR=/home/louism/Documents/VGVAPG
PYTHONDIR=/home/louism/Documents/VariationalApproach

# overlap matrix
python $PYTHONDIR/vp_correlationMatrix.py \
	-inDir	$WORKDIR/nobackup/PBVTrace/PBV_SD/rep1 $WORKDIR/nobackup/PBVTrace/PBV_SD/rep2 $WORKDIR/nobackup/PBVTrace/PBV_SD/rep3 $WORKDIR/nobackup/PBVTrace/PBV_SD/rep4 \
	-outDir $WORKDIR/correlationMatrix/PBV_SD \
	-PBVFile $WORKDIR/pbvsets/PBV_SDxxx.dat \
	-overlapMatrix \

# correlation matrix
for TAU in 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000; do
        python $PYTHONDIR/vp_correlationMatrix.py \
	-inDir	$WORKDIR/nobackup/PBVTrace/PBV_SD/rep1 $WORKDIR/nobackup/PBVTrace/PBV_SD/rep2 $WORKDIR/nobackup/PBVTrace/PBV_SD/rep3 $WORKDIR/nobackup/PBVTrace/PBV_SD/rep4 \
        	-outDir $WORKDIR/correlationMatrix/PBV_SD \
        	-PBVFile $WORKDIR/pbvsets/PBV_SDxxx.dat \
		-reversible \
		-lagTime $TAU
done

