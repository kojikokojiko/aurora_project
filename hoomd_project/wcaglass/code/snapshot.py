import sys
import os

import itertools
import math

import gsd.hoomd
import hoomd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches

import matplotlib.patches as pat
import sys
import os
import random
import time
import pickle
from PIL import Image

 

first_time=time.time()
# INPUT PARAMETER
nu=float(sys.argv[1])
fixed_percent=float(sys.argv[2])
kbT=float(sys.argv[3])
 
 

random.seed(12)
# this section is fixed  #############
#particle parameter
Nx = 64
Ny = 64
#粒子数
N = Nx*Ny 

lx = 64
ly = lx

#AL  半径比
AL = 1.4
EL = 1/3 #個数比
N1 = int(N*EL) #大粒子の個数　
N0 = N-N1 #小粒子の個数
##############################


m=1.0
epsilon=2.0


# 刻み幅小さすぎの可能性もあるから大きめにしてみてもいいかも

# 総時間幅
real_time=250
# 出力を始めるまでの時間幅＝＞平衡状態に達するまでの時間
after_time=20
# 出力をしている時間幅
pos_out_time=(real_time-after_time)
# 刻み時間
dt = 4.0e-5
# 総ステップ数
nsteps=int(real_time/dt)
# 出力を始めるまでのステップ数
after_steps=int(after_time/dt)
# 出力をしているときの送ステップ数
pos_out_steps=int(pos_out_time/dt)

########################################################出力間隔調整
# # 出力をしているときのステップ数の間隔
# pos_out_steps_period=5
# # 総出力回数
# num_out=int(pos_out_steps/pos_out_steps_period)
########################################################総出力回数調整
# 総出力回数
num_out=int(4000)
# 出力をしているときのステップ数の間隔
pos_out_steps_period=int(pos_out_steps/num_out)
##############################################################
# image出力枚数
image_out_num=2
# 総出力回数に対してイメージを出力させる間隔
image_out_period=int(num_out/image_out_num)


ver=str(nu)+"_"+str(fixed_percent)+"_"+str(kbT)

main_dir="./"+ver
if not os.path.exists(main_dir): os.makedirs(main_dir)
print(main_dir)



traj = gsd.hoomd.open(main_dir+"/log_pos.gsd",mode='r')






# traj = gsd.hoomd.open('log_force2d_'+ver+'.gsd', 'rb')


output_dir=main_dir+"/figure_2d"
if not os.path.exists(output_dir): os.makedirs(output_dir)
figsize=(10,10)
plt.figure(figsize=figsize,dpi=150)

print(len(traj))
sigma=traj[-1].particles.diameter
set_sigma=list(set(sigma))

typeid=traj[-1].particles.typeid
print(typeid)
print(sigma)
print(set_sigma)
for t in range(len(traj)-1,0,-image_out_period):
    print(t)
    bx=plt.axes()
    plt.axis([-lx/2,lx/2,-ly/2,ly/2])
    position=traj[t].particles.position
   
    for i in range(N):
        # print(i)

            # Circleパッチを作成
        if (sigma[i]==set_sigma[0]):
            c="r"
        else:
            c="b"

        c=pat.Circle(xy=(position[i][0],position[i][1]),radius=sigma[i]/2,fc=c)
        if i<100:
            plt.annotate(str(i),(position[i][0],position[i][1]))
        bx.add_patch(c)
    

    plt.title("step"+str(t))
    plt.savefig(output_dir+"/figure{0}.png".format(t))
    plt.cla()
    plt.clf()

    ############アニメーション################    
images=[]
# image_num=sum(os.path.isfile(os.path.join(pic_output name)) for name in os.listdir(pic_output))
image_num=sum(os.path.isfile(os.path.join(output_dir,name))for name in os.listdir(output_dir))
print(image_num)
for i in range(len(traj)-1,0,-image_out_period):
    file_name=output_dir+"/figure"+str(i)+".png"
    im=Image.open(file_name)
    images.append(im)

gif_output_dir=main_dir+"/abpgif2"

if not os.path.exists(gif_output_dir): os.makedirs(gif_output_dir)
images[0].save(gif_output_dir+"/out_ela2.gif",save_all=True,append_images=images[1:],loop=0,duration=10)
    
    