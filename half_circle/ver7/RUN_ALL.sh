#!/bin/bash

temp_rho_list=( 0.2  0.4  0.6  0.8)
v0_list=(0.1 1.0 2.0)
g++ -std=c++11 -o half_circle half_circle.cpp 
for rho in "${temp_rho_list[@]}"
do 
    for v0 in  "${v0_list[@]}"
    do
        qsub RUN $rho $v0
    done
done

# g++ -std=c++11 -o iabp2 iABP2.cpp
# nohup time ./iabp2 1.0 1.0 10000 # Pe M N
	






# rho=float(sys.argv[1])
# ave_flow=float(sys.argv[2])
# static_dia=float(sys.argv[3])
# reduced_speed=float(sys.argv[5])
# rotational_diffusion=float(sys.argv[6])

