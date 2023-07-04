import sys
import os
sys.path.append('/home/isobelab2022/build3/hoomd')
import itertools
import math
from scipy.spatial import Voronoi, voronoi_plot_2d
from shapely.geometry import Polygon
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
Lx=L
Ly=L

media_dir="/media/isobelab2022/data/normal_glass/ver2"



ver=str(nu)+"_"+str(fixed_percent)+"_"+str(kbT)
main_dir="./"+ver

vornoi_dir=main_dir+"/vornoi"
if not os.path.exists(vornoi_dir):
    os.makedirs(vornoi_dir)

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

shifts = np.array([[0, 0], [0, Ly], [0, -Ly], [Lx, 0], [-Lx, 0], [Lx, Ly], [Lx, -Ly], [-Lx, Ly], [-Lx, -Ly]])



temp_lx=35
temp_ly=35


plt.figure(figsize=(10,10))


def calculate_area(vor, index):
    """Calculate the area of a Voronoi cell by index."""
    # ボロノイセルの頂点を取得
    vertices = vor.regions[vor.point_region[index]]
    if -1 in vertices:  # 無限遠の点を含むセルは無視
        return np.inf
    else:
        # 頂点を使用して多角形を作成し、その面積を計算
        return Polygon([vor.vertices[i] for i in vertices]).area


vornoi_list=[]
for t in range(T):
    points=all_pos[t]
    
    points_periodic = np.vstack([points + shift for shift in shifts])

    # 条件に合致するものだけを抽出
    mask = (points_periodic[:, 0] >= -temp_lx) & (points_periodic[:, 0] <= temp_lx) & (points_periodic[:, 1] >= -temp_ly) & (points_periodic[:, 1] <= temp_ly)
    points_filtered = points_periodic[mask]
    # ボロノイ図を計算
    vor = Voronoi(points_filtered)
    

    # ボロノイ図をプロット
    voronoi_plot_2d(vor)
    # 原点を中心にプロットを制限
    plt.xlim(-Lx/2, Lx/2)
    plt.ylim(-Ly/2, Ly/2)
    
    plt.savefig(main_dir+"/voronoi_"+str(t)+".png")
    
    plt.clf()
    # 各セルの面積を計算
    areas = [calculate_area(vor, i) for i in range(len(points_filtered))]

    vornoi_list.append(areas)
    

vornoi_list=np.array(vornoi_list)

    
    
            
            
    

np.savez(main_dir+"/voronoi.npz",chi4_t=vornoi_list,zerod_time=zerod_time)

print("FiNISH")

del traj