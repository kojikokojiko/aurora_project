#!/bin/bash



rho_list=(0.8)
ave_flow_list=(1.0 3.0 5.0 )
static_dia=50.0
reduced_speed_list=(10.0 )
rotational_diffusion_list=(0.1 10.0 )
for rho in "${rho_list[@]}"
do 
    for ave_flow in  "${ave_flow_list[@]}"
    do
        for reduced_speed in  "${reduced_speed_list[@]}"
        do

            for rotational_diffusion in  "${rotational_diffusion_list[@]}"
            do
                qsub RUN $rho $ave_flow $static_dia $reduced_speed  $rotational_diffusion
            done
            
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

