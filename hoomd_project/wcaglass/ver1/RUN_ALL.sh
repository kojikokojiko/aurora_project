#!/bin/bash



nu_list=(0.78  )
# 0.55 0.6 0.65 0.7 0.75 0.8)

fixed_per_list=(0.0 0.1  0.3  0.5)
kbT_list=(0.5 1.0 1.5)
for nu in "${nu_list[@]}"
do 
    for kbT in  "${kbT_list[@]}"
    do
        for fixed_per in  "${fixed_per_list[@]}"
        do
        echo $nu $kbT  $fixed_per
        name="${nu}_${kbT}_${fixed_per}"
        qsub RUN $nu $fixed_per $kbT
        # qsub RUN_SNAP $nu $fixed_per $kbT
        # nohup time python animation.py $nu $pl $fixed_per> ${name}.log 2>&1 &
        done
    done
done

# g++ -std=c++11 -o iabp2 iABP2.cpp
# nohup time ./iabp2 1.0 1.0 10000 # Pe M N
	






# rho=float(sys.argv[1])
# ave_flow=float(sys.argv[2])
# static_dia=float(sys.argv[3])
# reduced_speed=float(sys.argv[5])
# rotational_diffusion=float(sys.argv[6])

