#!/bin/bash

for i in Dis*; do
	echo $i;
	cat $i/Traj_info_Dist*.xyz >> train.xyz;
done
