#!/bin/sh

#for RES in A C D E F G H I K L M N P Q R S T V W Y; do
#	mkdir Ac_$RES_NHMe	
#	#mv Ac\_$RES\_NHMe.dat Ac\_$RES\_NHMe/rEV\_norm.dat
#done

for RES in A C D E F G H I K L M N P Q R S T V W Y; do
	mkdir Ac\_$RES\_NHMe
	cp Ac\_$RES\_NHMe.dat Ac\_$RES\_NHMe/rEV\_norm.dat
done
