#!/bin/sh

WORKDIR=/home/louism/Documents/Ac_A_NHMe
PYTHONDIR=/home/louism/Documents/VariationalApproach

# overlap matrix
for TAU in  000010 000020 000030 000040 000050 000060 000070 000080 000090 000100 000200 000300 000400 000500 000600 000700 000800 000900 001000; do
	python $PYTHONDIR/vp_generalizedEV.py \
		-S $WORKDIR/correlationMatrix/PBVs_S/S.npy	\
		-C $WORKDIR/correlationMatrix/PBVs_S/C_tau$TAU.npy \
		-eigenValueFile $WORKDIR/correlationMatrix/PBVs_S/l_tau$TAU.dat \
		-eigenVectorFile $WORKDIR/correlationMatrix/PBVs_S/a_tau$TAU.dat \
		-ascii 
done
