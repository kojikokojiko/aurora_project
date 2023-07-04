import sys
import os
sys.path.append('/home/isobelab2022/build3/hoomd')
import itertools
import math

import gsd.hoomd
import hoomd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from numba import jit, f8, i8, b1, void,njit
import matplotlib.patches as pat
import sys
import os
import fresnel
from PIL import Image
import IPython
import packaging.version
import random
import time
import pickle


nu=float(sys.argv[1])
fixed_percent=float(sys.argv[2])
kbT=float(sys.argv[3])
L=64
LX=L
LY=L

media_dir="/media/isobelab2022/data/normal_glass/ver2"



ver=str(nu)+"_"+str(fixed_percent)+"_"+str(kbT)
main_dir="./"+ver


traj_dir=media_dir+"/"+ver
traj_path=traj_dir+"/log_pos_"+ver+".gsd"
traj = gsd.hoomd.open(traj_path, 'rb')

data_path=media_dir+"/"+ver+"/data.pickle"
with open(data_path, mode='rb') as f:
    data = pickle.load(f)



dt=data["dt"]
pos_out_steps_period=data["pos_out_steps_period"]
small_sigma=data["sigma_ss"]
msd_times=[i*dt*pos_out_steps_period for i in range(len(traj)-1)]
zerod_time=msd_times/small_sigma

zerod_time=np.array(zerod_time)





all_pos=[traj[i].particles.position for i in range(len(traj)-1)]


all_pos=np.array(all_pos)
# posは（時間、粒子数、3次元[x,y,z]）の次元になっている。


N = all_pos.shape[1]
T=all_pos.shape[0]

print(N,T)




def w(distance):
    if distance<=0.3*small_sigma:
        return distance
    else:
        return 0




def periodic_distance(r1, r2, L):
    # Calculate the distance between two points in a periodic boundary system
    dr = r2 - r1
    dr = dr - L * np.round(dr / L)
    return np.sqrt(np.sum(dr**2))



def Q(t, r, L):
    N = len(r)
    sum_ = 0
    for j in range(N):
        for i in range(N):
            distance = periodic_distance(r[j, 0], r[i, t], L)
            sum_ += w(distance)
    return sum_/N




Qt=[]
chi4_t=[]
for t in range(T):
    Q_t =Q(t,all_pos[t],L)
    Q_2_t =Q(t,all_pos[t],L)**2
    
    Q_2_t_avg = np.mean(Q_2_t) # This needs to be calculated from multiple realizations
    Q_t_avg = np.mean(Q_t) # This needs to be calculated from multiple realizations
    
    chi4 = (V /(N**2 *kbT))* (Q_2_t_avg - Q_t_avg**2)
    
    chi4_t.append(chi4)
    
    
            
            
    
chi4_t=np.array(chi4_t)
np.savez(main_dir+"/chi4_t.npz",chi4_t=chi4_t,zerod_time=zerod_time)

print("FiNISH")

del traj