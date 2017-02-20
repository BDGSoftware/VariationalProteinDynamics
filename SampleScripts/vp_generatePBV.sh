#!/bin/sh

OUTDIR=/home/louism/Documents/Ac_AV_NHMe/pbvsets
PYTHONDIR=/home/louism/Documents/VariationalApproach

# SD
python $PYTHONDIR/vp_PBVBasisSet.py \
-outDir $OUTDIR \
-sequence $OUTDIR/sequence.dat \
-single \
-double 
