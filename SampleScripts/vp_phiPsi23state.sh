#!/bin/sh

WORKDIR=/home/louism/Documents/Ac_A_NHMe
PYTHONDIR=/home/louism/Documents/VariationalApproach

for REP in rep1 rep2 rep3 rep4 rep5 rep6 rep7 rep8 rep9 rep10 rep11 rep12 rep13 rep14 rep15 rep16 rep17 rep18 rep19 rep20; do
	python  $PYTHONDIR/vp_phiPsi23state.py \
		-inDir  $WORKDIR/phiPsiTrace/$REP \
		-outDir $WORKDIR/nobackup/microstateTrace/PBVs/$REP \
		-sequence $WORKDIR/pbvsets/sequence.dat \
		-shift  \
		-tolerance 0.0001 
done
