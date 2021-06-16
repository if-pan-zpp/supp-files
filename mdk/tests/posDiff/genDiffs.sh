#!/bin/bash

for t in 1ubq 9aac
do
	echo $t
	cp cg.f $t/data
	cd $t/data
	patch cg.f cg_patch
	echo "Compiling..."
	gfortran -O3 cg.f -o cg 2> /dev/null
	echo "Running..."
	./cg inputfile > positions.txt
	rm discard*
	cd ../..
done
