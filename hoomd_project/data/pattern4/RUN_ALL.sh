#!/bin/bash



rho_list=(0.8)
ave_flow_list=( 0.0)
static_dia=150.0
reduced_speed_list=(0.0 )
rotational_diffusion_list=(5.0)
for rho in "${rho_list[@]}"
do 
    for ave_flow in  "${ave_flow_list[@]}"
    do
        for reduced_speed in  "${reduced_speed_list[@]}"
        do

            for rotational_diffusion in  "${rotational_diffusion_list[@]}"
            do
                echo $rho $ave_flow $static_dia $reduced_speed  $rotational_diffusion
             
                qsub RUN $rho $ave_flow $static_dia $reduced_speed  $rotational_diffusion 
           
            done
            
        done
    done
done
