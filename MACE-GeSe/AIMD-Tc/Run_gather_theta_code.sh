#!/bin/bash

for i in {200..700..25}
do
    ( echo $i
      cd Test_Temp$i
      echo pwd
      python Gather_Data_from_traj-updated.py > Theta2.out
    )  
done
