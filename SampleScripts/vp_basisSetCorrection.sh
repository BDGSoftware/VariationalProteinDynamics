#!/bin/sh

WORKDIR=/home/louism/Documents/hIAPP
PYTHONDIR=/home/louism/Documents/VariationalApproach/variational

python -i $PYTHONDIR/vp_basisSetCorrection.py \
	-inDir	$WORKDIR/correlationMatrix/PBV_SD_full \
	-outDir $WORKDIR/correlationMatrix/PBV_SDTest \
	-PBVFile $WORKDIR/pbvsets/PBV_SDxxx_full.dat \
	-PBVFileOut $WORKDIR/pbvsets/PBV_SDxxxC.dat
