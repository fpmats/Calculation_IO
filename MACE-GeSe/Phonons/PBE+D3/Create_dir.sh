#!/bin/bash

for i in {01..08}
do
	mkdir Force-$i
	cp POSCAR-0$i Force-$i/POSCAR
	cp ASE-VASP-SCF.py Force-$i/.
	cp run_VASP-hprthread.sh Force-$i/.
done
