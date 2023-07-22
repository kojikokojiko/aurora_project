import sys
import os

import matplotlib.pyplot as plt
import matplotlib.patches as pat
import gsd.hoomd
import sys
from PIL import Image
import os
import math
import  numpy as np
import sys





rho=float(sys.argv[1])
ave_flow=float(sys.argv[2])
static_dia=float(sys.argv[3])
reduced_speed=float(sys.argv[4])
rotational_diffusion=float(sys.argv[5])
ver=str(rho)+"_"+str(ave_flow)+"_"+str(static_dia)+"_"+str(reduced_speed)+"_"+str(rotational_diffusion)



main_dir="./"+ver


# temp_dir="./"+ver+"/log_pos_"+ver+".gsd"
temp_dir="./"+ver+"/log_pos.gsd"
phi6_2_dir=main_dir+"/phi6_2_scatter"
# i_phi_dir=main_dir+"/i_phi"
if not os.path.exists(phi6_2_dir): os.makedirs(phi6_2_dir)
# if not os.path.exists(i_phi_dir): os.makedirs(i_phi_dir)


print(temp_dir)
traj = gsd.hoomd.open(temp_dir)
# traj = gsd.hoomd.open(dir, 'rb')




lx=traj[0].configuration.box[0]
ly=traj[0].configuration.box[1]
ppi=72 # points per inche 

plt.figure(figsize=(10,10*ly/lx))

NP=len(traj[0].particles.position)
print(NP)

sigma=1.0
# あんま頭よくない処理
size=sigma


output_list=[]
phi6_output_list=[]
# 描画設定##############################################
# https://qiita.com/stanaka2/items/c40841f858d7083aad4e
fig = plt.figure(figsize=(10,10*ly/lx), dpi=300.0)
ax  = fig.add_axes((0.1,  0.1, 0.8, 0.8))

# 枠の範囲指定
xmin,xmax=-lx/2.0, lx/2.0
ymin,ymax=-ly/2.0, ly/2.0

# 描写範囲の長さを取得(dpi単位)
# x軸をベースに計算しているがy軸でも同じ。アスペクト比が違えばおかしくなる
ax_length=ax.bbox.get_points()[1][1]-ax.bbox.get_points()[0][1]

# dpi単位の長さをポイント単位に変換(dpiからインチ単位にし、インチからポイント単位に変換)
ax_point = ax_length*ppi/fig.dpi

# x軸の実スケールをポイント単位に変換
xsize=xmax-xmin
fact=ax_point/xsize

# scatterのマーカーサイズは直径のポイントの二乗を描くため、実スケールの半径をポイントに変換し直径にしておく
size*=2*fact

####################################################

data=np.load(main_dir+"/phi6_2.npy")
data_index=0
# ここはnpyに出力されたforに合わせる必要がある。  
for t in range(len(traj)-1,len(traj)-10,-2):
    phi6_2=data[data_index]
    data_index+=1

    
    print(t)
    pos=traj[t].particles.position.T
    rx=pos[0]
    ry=pos[1]


    # extract_index1=np.where(phi6_2>0.0)
    # extract_phi61=phi6_2[extract_index1]
    # extract_rx1=rx[extract_index1]
    # extract_ry1=ry[extract_index1]


    # extract_index2=np.where(0.2<phi6_2)
    # extract_phi62=phi6_2[extract_index2]
    # extract_rx2=rx[extract_index2]
    # extract_ry2=ry[extract_index2]


    # extract_index3=np.where(0.55<phi6_2)
    # extract_phi63=phi6_2[extract_index3]
    # extract_rx3=rx[extract_index3]
    # extract_ry3=ry[extract_index3]

    # extract_index4=np.where(0.6<phi6_2)
    # extract_phi64=phi6_2[extract_index4]
    # extract_rx4=rx[extract_index4]
    # extract_ry4=ry[extract_index4]


    extract_index5=np.where(0.75<phi6_2)
    extract_phi65=phi6_2[extract_index5]
    extract_rx5=rx[extract_index5]
    extract_ry5=ry[extract_index5]
    
# 描画############################
# Real_phi6#################################
# dpiは任意の値を設定できる。デフォルトは100
    fig = plt.figure(figsize=(12.5,5.0), dpi=300.0)
    ax  = fig.add_axes((0.1,  0.1, 0.8, 0.8))

    ax.set_aspect("equal") # アスペクト比を等しくする
                        # 今は枠のサイズ設定で最初から等しいのであってもなくてもよい

    # 二乗にして与える 
    im=ax.scatter(rx,ry,s=size**2,c=phi6_2,cmap="cool", linewidths=0)
    im5=ax.scatter(extract_rx5,extract_ry5,s=size**2,c="red", linewidths=0)


    # for i in range(0,NP,50):
    #     c=pat.Circle(xy=(rx[i],ry[i]),radius=r_sann_1[i],fc="g",alpha=0.2)
    #     ax.add_patch(c)

    # for i in range(100,200):
    #     ax.text(rx[i], ry[i], i,fontsize="small")

    c=pat.Circle(xy=(lx/4-lx/2,ly/2-ly/2),radius=static_dia/2,fc="b")
    ax.add_patch(c) 
    # plt.colorbar(im,ticks=[0.0,0.2,0.4 ,0.6,0.75,0.8])

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    ax.set_xticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 5))
    ax.set_yticks(np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 5))

    # ax.grid(which='both', axis='both')

    # plt.colorbar()
    # plt.colorbar(sc, ax=ax)
    plt.title("phi6_2")
    plt.savefig(phi6_2_dir+"/phi6_2_{0}.png".format(t))
    # plt.savefig(main_dir+"/figure.png")
    # plt.cla()
    plt.close()
