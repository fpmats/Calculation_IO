#!/bin/bash

for i in {0..700..50}
do
    echo $i
    mkdir Test_Temp$i
    cp POSCAR Test_Temp$i/POSCAR
    cp md_stage.py Test_Temp$i/.
    cp model_swa.model Test_Temp$i/.
    cp job.sh Test_Temp$i/.
    #cp Gather* Test_Temp$i/.
    sed -i 's/temperature_K=400/temperature_K='$i'/g' Test_Temp$i/md_stage.py
done
