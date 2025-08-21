#!/bin/bash

#phonopy -d --dim 6 6 2 --pa auto

#phonopy -f Force-{01..08}/vasprun.xml

#phonopy-load phonopy_disp.yaml --mesh 31 31 31 -t -p -s

#phonopy-load phonopy_disp.yaml --band "0.0 0.0 0.0  0.5 0.0 0.0  0.5 0.5 0.0  0.0 0.5 0.0  0.0 0.0 0.0  0.0 0.0 0.5  0.5 0.0 0.5  0.5 0.5 0.5  0.0 0.5 0.5 0.0 0.0 0.5" --band-labels "$\Gamma$ X S Y $\Gamma$ Z U R T Z"  --mesh 31 31 31 --pdos "1, 2" -p -s
phonopy-load phonopy_disp.yaml --pa="1 0 0 0 1 0 0 0 1" --band "0.0 0.0 0.0  0.5 0.0 0.0  0.5 0.5 0.0  0.0 0.5 0.0  0.0 0.0 0.0  0.0 0.0 0.5  0.5 0.0 0.5  0.5 0.5 0.5  0.0 0.5 0.5 0.0 0.0 0.5" --band-labels "$\Gamma$ X S Y $\Gamma$ Z U R T Z"  --mesh 31 31 31 --pdos "1, 2" -p -s
