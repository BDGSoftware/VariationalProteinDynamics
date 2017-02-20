#!/bin/sh

WORKDIR=/home/louism/Documents/HSSNNFGAIL
PYTHONDIR=/home/louism/Documents/VariationalApproach/variational

python $PYTHONDIR/vp_basisVectors.py \
                -inDir  $WORKDIR/nobackup/microstateTrace/PBV \
                -outDir $WORKDIR/pbvsets/basisVectors/PBV \
                -sequence $WORKDIR/pbvsets/sequence.dat \
		-basisSet $PYTHONDIR/BasisVectors
