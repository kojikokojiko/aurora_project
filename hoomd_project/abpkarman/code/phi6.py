import sys
import os

import matplotlib.pyplot as plt
import matplotlib.patches as pat
import gsd.hoomd
import hoomd
import sys
from PIL import Image
import os
import math
import  numpy as np
import seaborn as sns
import sys


def make_neighbor_list(around_grid_num):
    if around_grid_num==5:
        neighbor_col=[
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1,    1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
            -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
        ]
        neighbor_row=[
            -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,
            -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
            -3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,
            -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
            -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
             0, 0, 0, 0, 0,    0, 0, 0, 0, 0,
             1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
             2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
             3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
             5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5

        ]
    return neighbor_row,neighbor_col





rho=float(sys.argv[1])
ave_flow=float(sys.argv[2])
static_dia=float(sys.argv[3])
reduced_speed=float(sys.argv[4])
rotational_diffusion=float(sys.argv[5])
ver=str(rho)+"_"+str(ave_flow)+"_"+str(static_dia)+"_"+str(reduced_speed)+"_"+str(rotational_diffusion)



main_dir="./"+ver




phi6_2_dir=main_dir+"/phi6_2"
# i_phi_dir=main_dir+"/i_phi"
if not os.path.exists(phi6_2_dir): os.makedirs(phi6_2_dir)
# if not os.path.exists(i_phi_dir): os.makedirs(i_phi_dir)


print(main_dir)
traj = gsd.hoomd.open(main_dir+"/log_pos.gsd")
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

around_grid_num=5
neighbor_row,neighbor_col= make_neighbor_list(around_grid_num)
neighbor_len=len(neighbor_row)


bo=(lx/0.4)
l_gx_d=lx/bo
l_gx=l_gx_d
l_gy=l_gx
n_gx=math.ceil(lx/l_gx)
n_gy=math.ceil(ly/l_gy)
n_g_len=n_gx*n_gy

r_cut_sann=l_gx_d*around_grid_num
pair_length_g=50
# NN_MAX=20
NN_MAX=40
Nmax=80


# 検査体積 binn
l_grid_bin=1.5
n_grid_bin=math.ceil(lx/l_grid_bin)
Vc=l_grid_bin*ly
bin_x=[i*l_grid_bin for i in range(n_grid_bin)]


output_list=[]
phi6_output_list=[]
# 描画設定##############################################
# https://qiita.com/stanaka2/items/c40841f858d7083aad4e
fig = plt.figure(figsize=(12.5,5.0), dpi=100.0)
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
for t in range(len(traj)-1,len(traj)-20,-2):
    
    print(t)
    pos=traj[t].particles.position.T
    rx=pos[0]
    ry=pos[1]


    # Gmap_create########################
    g_map=np.full((n_gy,n_gx),-1,dtype="int")
    pair_list_g=np.full((NP,pair_length_g),-1,dtype="int")
    for i in range(NP):
        gx_map=int((rx[i]+lx/2)/l_gx)
        gy_map=int((ry[i]+ly/2)/l_gy)
        if(g_map[gy_map][gx_map]!=-1):
            print("TOO LARGE GRID")
            sys.exit()
        g_map[gy_map][gx_map]=i
    # np.savetxt(main_dir+"/g_map{0}.txt".format(t),g_map)

    ##################################################################
    # make_pairlist_g####################################

    for i in range(n_gy):
        for j in range(n_gx):
            select_index=g_map[i][j]
            if select_index==-1:
                continue
            particle_counter=0
            for k in range(neighbor_len):
                search_gx=j+neighbor_col[k]
                search_gy=i+neighbor_row[k]
                if(search_gx>=n_gx):
                    search_gx-=n_gx
                elif(search_gx<0):
                    search_gx+=n_gx
                if (search_gy>=n_gy):
                    search_gy-=n_gy
                elif(search_gy<0):
                    search_gy+=n_gy
                
                search_index=g_map[search_gy][search_gx]
                if search_index==select_index:
                    print("what wrong")
                    print(select_index)
                if search_index!=-1:
                    pair_list_g[select_index][particle_counter]=search_index
                    particle_counter+=1

            if particle_counter==pair_length_g-1:
                print("ERROR HAPPEN :PAIR_LENGH_G is few")
                sys.exit()
                
            pair_list_g[select_index][-1]=particle_counter
    ############################################
 

    # np.savetxt(main_dir+"/pair_list{0}.txt".format(t),pair_list_g)


    # SANN_CUT###########################
    NNnum=np.full(NP,0,dtype="int")
    NN=np.zeros((NP,NN_MAX),dtype="int")
    NNnum_S=np.zeros(NP,dtype="int")
    NN_S=np.zeros((NP,NN_MAX),dtype="int")
    DR_S=np.zeros((NP,NN_MAX))
    for i in range(NP):
        pair_count=0
        loop_length=pair_list_g[i][-1]

        for loop_count in range(loop_length):
            pair_index=pair_list_g[i][loop_count]
            rxij=rx[i]-rx[pair_index]
            ryij=ry[i]-ry[pair_index]
            if(rxij>=lx/2):
                rxij=rxij-lx
            elif (rxij<=-lx/2):
                rxij=rxij+lx

            if(ryij>=ly/2):
                ryij=ryij-ly
            elif (ryij<=-ly/2):
                ryij=ryij+ly
            r2_sann=rxij*rxij+ryij*ryij
            dr_sann=np.sqrt(r2_sann)
            if(dr_sann<=r_cut_sann):
                NN_S[i][NNnum_S[i]]=pair_index
                DR_S[i][NNnum_S[i]]=dr_sann
                NNnum_S[i]+=1



    ########################################################


    # ???????????????????なんで―1????????????????????????????????????????
    # for i in range(N):
    #     NNnum_S[i]=NNnum_S[i]-1
    ################################

    # SANNN SORT###############################################
    for i in range(NP):
        NNs_list=[[DR_S[i][k],NN_S[i][k]] for k in range(NNnum_S[i])]
        # ???????????????????なんで＋１
        # NNs_list=[[DR_S[i][k],NN_S[i][k]] for k in range(NNnum_S[i]+1)]
        NNs_list.sort(key=lambda x:x[0])
        for k in range(NNnum_S[i]):
            DR_S[i][k]=NNs_list[k][0]
            NN_S[i][k]=NNs_list[k][1]

    # np.savetxt(main_dir+"/NNnum_S{0}.txt".format(t),NNnum_S)
    # np.savetxt(main_dir+"/NN_S{0}.txt".format(t),NN_S)
    # np.savetxt(main_dir+"/DR_S{0}.txt".format(t),DR_S)

    # MAIN_SANN#########################################################
    r_sann_1=np.zeros(NP)
    for i in range(NP):
        for j in range(3,NNnum_S[i]):
            sann_m=j
            sr_num=NNnum_S[i]
            r_min=DR_S[i][sann_m-1]
            r_max=DR_S[i][sr_num-1]
            for n in range(Nmax):
             
                r_mid=(r_min+r_max)/2
                sum_rm_min=0.0
                sum_rm_max=0.0
                sum_rm_mid=0.0
                for m in range(sann_m):
                    sum_rm_min=sum_rm_min+math.acos(DR_S[i][m]/r_min)
                    sum_rm_max=sum_rm_max+math.acos(DR_S[i][m]/r_max)
                    sum_rm_mid=sum_rm_mid+math.acos(DR_S[i][m]/r_mid)
                f_max=sum_rm_max-np.pi
                f_min=sum_rm_min-np.pi
                f_mid=sum_rm_mid-np.pi
                if abs(f_mid)<1.0e-10:
                    break
                if f_mid*f_min>0.0:
                    r_min=r_mid
                    f_min=f_mid
                else:
                    r_max=r_mid
                    f_max=f_mid
            r_sann_m=r_mid
            if r_sann_m <=DR_S[i][sann_m]:
                NNnum[i]=sann_m
                r_sann_1[i]=r_sann_m
                for k in range(NNnum[i]):
                    NN[i][k]=NN_S[i][k]
                break
            else:
                sann_m=sann_m+1

    # np.savetxt(main_dir+"/NNnum{0}.txt".format(t),NNnum)
    # np.savetxt(main_dir+"/NN{0}.txt".format(t),NN)
    # np.savetxt("./DR.txt",DR_S)

    ###################################################################
    r_phi6=np.zeros(NP)
    i_phi6=np.zeros(NP)
    for i in range(NP):
        for j in range(NNnum[i]):
            NN_index=NN[i][j]
            dx=rx[NN_index]-rx[i]
            dy=ry[NN_index]-ry[i]

            if(dx>=lx/2):
                dx=dx-lx
            elif(dx<-lx/2):
                dx=dx+lx
            if(dy>=ly/2):
                dy=dy-ly
            elif(dy<-ly/2):
                dy=dy+ly
            r2_NN=dx*dx+dy*dy
            r_NN=math.sqrt(r2_NN)
            # if(r_NN==0.0):
                # print("r_NN is 0")
                # print(t)
                # print(i)
                # print(NN_index)
            temp=dx/r_NN
            if(temp>1.0):
                # print("DANGER")
                theta=0
            elif (temp<-1.0):
                # print("DANGER")/
                theta=math.pi
            else:
                if(dy>=0.0):
                    theta=math.acos(temp)
                elif(dy<0):
                    theta=-math.acos(temp)
            r_phi6[i]+=math.cos(6*theta)
            i_phi6[i]+=math.sin(6*theta)

    for i in range(NP):
        Ni=NNnum[i]
        if Ni==0:
            # print(r_phi6[i])
            # print(i_phi6[i])
            continue
            
        r_phi6[i]/=Ni
        i_phi6[i]/=Ni


    phi6_2=np.sqrt(r_phi6*r_phi6+i_phi6*i_phi6)

    phi6_output_list.append(phi6_2)







    phi6_2_list_bin=np.zeros(n_grid_bin)
    partice_count_list_bin=np.zeros(n_grid_bin)
    for temp_index in range(NP):
        gx_id=int((rx[temp_index]+lx/2.0)/l_grid_bin)
        phi6_2_list_bin[gx_id]+=phi6_2[temp_index]
        partice_count_list_bin[gx_id]+=1





    # print(partice_count_list_bin)
    for temp_gx_id in range(n_grid_bin):
        if partice_count_list_bin[temp_gx_id]==0:
            print("WARNING")
            print(temp_gx_id)
        phi6_2_list_bin[temp_gx_id]/=partice_count_list_bin[temp_gx_id]

    output_list.append(phi6_2_list_bin)
    print(len(output_list))


# 描画############################
# Real_phi6#################################
# dpiは任意の値を設定できる。デフォルトは100
    fig = plt.figure(figsize=(12.5,5.0), dpi=100.0)
    ax  = fig.add_axes((0.1,  0.1, 0.8, 0.8))

    ax.set_aspect("equal") # アスペクト比を等しくする
                        # 今は枠のサイズ設定で最初から等しいのであってもなくてもよい

    # 二乗にして与える 
    im=ax.scatter(rx,ry,s=size**2,c=phi6_2,cmap="cool", linewidths=0)
    # for i in range(0,NP,50):
    #     c=pat.Circle(xy=(rx[i],ry[i]),radius=r_sann_1[i],fc="g",alpha=0.2)
    #     ax.add_patch(c)

    # for i in range(100,200):
    #     ax.text(rx[i], ry[i], i,fontsize="small")

    c=pat.Circle(xy=(lx/4-lx/2,ly/2-ly/2),radius=static_dia/2,fc="b")
    ax.add_patch(c) 
    plt.colorbar(im)

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





# # Real_phi6#################################
#     plt.scatter(rx,ry,s=size**2,c=r_phi6,cmap="seismic", linewidths=0,vmin=-1.0,vmax=1.0)
#     plt.colorbar()
#     plt.savefig(r_phi_dir+"/r_phi_color{0}.png".format(t))
#     # plt.clf()
#     # plt.cla()
#     plt.close()


output_list=np.array(output_list)
np.save(main_dir+"/phi6_bin2.npy", output_list)

phi6_output_list=np.array(phi6_output_list)
np.save(main_dir+"/phi6_2.npy", phi6_output_list)






# # Imagenary_phi6#################################
# # dpiは任意の値を設定できる。デフォルトは100
#     fig = plt.figure(figsize=(12.5,5.0), dpi=100.0)
#     ax  = fig.add_axes((0.1,  0.1, 0.8, 0.8))

#     ax.set_aspect("equal") # アスペクト比を等しくする
#                         # 今は枠のサイズ設定で最初から等しいのであってもなくてもよい

#         # 二乗にして与える 
#     ax.scatter(rx,ry,s=size**2,c=i_phi6,cmap="seismic", linewidths=0,vmin=-1.0,vmax=1.0)
#     # plt.colorbar()

#     ax.set_xlim(xmin,xmax)
#     ax.set_ylim(ymin,ymax)

#     ax.set_xticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 5))
#     ax.set_yticks(np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 5))

#     # ax.grid(which='both', axis='both')

#     # plt.colorbar()
#     # plt.colorbar(sc, ax=ax)
#     plt.title("i_phi6")
#     plt.savefig(i_phi_dir+"/i_phi{0}.png".format(t))

#     # plt.cla()
#     plt.close()

# # get cmap############################################################
#     plt.scatter(rx,ry,s=size**2,c=i_phi6,cmap="seismic", linewidths=0,vmin=-1.0,vmax=1.0)
#     plt.colorbar()
#     plt.savefig(i_phi_dir+"/i_phi_color{0}.png".format(t))
#     # plt.cla()
#     plt.close()


####################################################################


    # local_phi6=0.0
    # for i in range(N):
    #     local_phi6+=math.sqrt(r_phi6[i]**2+i_phi6[i]**2)
    # local_phi6/=N

    # mean_re=np.mean(r_phi6)
    # mean_im=np.mean(i_phi6)
    # global_phi6=math.sqrt(mean_re**2+mean_im**2)

    # global_phi6_list[t]=global_phi6
    # local_phi6_list[t]=local_phi6



# def sun_cut(NP,r_cut_sann,rx,ry,lx,ly,pair_list_g):
#     pair_length_sann=50
#     dr_list_sann=np.zeros(NP,pair_length_sann)
#     pair_index_sann=np.full((NP,pair_length_sann),-1)
#     for i in range(NP):
#         pair_count=0
#         loop_length=pair_list_g[i][-1]

#         for loop_count in range(loop_length):
#             pair_index=G_MAP[i][loop_count]
#             rxji=rx[j]-rx[i]
#             ryji=ry[j]-ry[i]
#             if(rxji>=lx/2):
#                 rxji=rxji-lx
#             elif (rxji<=-lx/2):
#                 rxji=rxji+lx

#             if(ryji>=ly/2):
#                 ryji=ryji-ly
#             elif (ryji<=-ly/2):
#                 rxji=ryij+ly
#             r2=rxji*rxji+ryij*ryij
#             if(r2<=r_cut_sann**2):
#                 pair_index_sann[i][k]=pair_index
#                 dr_list_sann[i][k]=np.sqrt(r2)
#                 pair_count+=1
#         dr_list_sann[i][-1]=pair_count

            
#     return dr_list_sann,pair_index_sann
    
# def sann(dr_list_sann,pair_index_sann,NP):
#     N_MAX=50
#     NN_NUM=np.zeros(N_MAX)
#     NN=np.zeros(NP,N_MAX)

#     for i in range(NP):
#         num_pair=pair_index_sann[i][-1]
#         for j in range(3,num_pair):
#             SANN_M=j
#             r_min=dr_list_sann[i][SANN_M]
#             r_max=dr_list_sann[i][num_pair]
#             for n in range(N_MAX):
#                 r_mid=(r_min+r_max)/2
#                 sum_rm_min=0.0
#                 sum_rm_max=0.0
#                 sum_rm_mid=0.0
#                 for m in range(SANN_M):
#                     sum_rm_min+=math.acos(dr_list_sann[i][m]/r_min)
#                     sum_rm_max+=math.acos(dr_list_sann[i][m]/r_max)
#                     sum_rm_mid+=math.acos(dr_list_sann[i][m]/r_mid)
#                 f_min=sum_rm_min-math.pi
#                 f_max=sum_rm_max-math.pi
#                 f_mid=sum_rm_max-math.pi
#                 if(abs(f_mid)<=1e-10):
#                     r_sann_m=dr_mid
#                     if(r_sann_m<=dr_list_sann[i][SANN_M+1]):
#                         NN_NUM[i]=SANN_M
#                         for k in range(NN_NUM[i]):
#                             NN[i][k]=NN_S[i][k]
                        
#                 if(f_mid*f_min>=0.0):
#                     r_min=r_mid
#                     f_min=f_mid
#                 else:
#                     r_max=r_mid
#                     f_max=f_mid







# def sann_sort(r_tmp,dr,NN_S,NN_NUM_S):
#     for i in range(N):
#         sr_num=int(NN_NUM_S[i])
#         for j in range(NN_MAX):
#             r_tmp=dr[i][j]
#             i_tmp=NN_S[i][j]
#         im_inssrot()
#         for j in range(NN_MAX):
#             dr[i][j]=R_





# def make_neighbor_list(around_grid_num):
#     if around_grid_num==5:
#         neighbor_col=[
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1,    1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
#         ]
#         neighbor_row=[
#             -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,
#             -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
#             -3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,
#             -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
#             -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
#              0, 0, 0, 0, 0,    0, 0, 0, 0, 0,
#              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#              2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
#              3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
#              4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
#              5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5

#         ]


#     # neighbor_row=[]
#     # neighbor_row_value=-around_grid_num
#     # for i in range(around_grid_num+1):
#     #     for j in range(around_grid_num):
#     #         neighbor_row.append(neighbor_row_value)
#     #     neighbor_row_value+=1
    
#     # for i in range(around_grid_num):
#     #     for j in  range(around_grid_num+1):
#     #         neighbor_row.append(neighbor_row_value)
#     #     neighbor_row_value+=1
    
#     # # yの近接リスト作成
#     # neighbor_col=[]
#     # neighbor_col_value=-around_grid_num
#     # for i in range(0,around_grid_num+1):
#     #     for j in range(around_grid_num):
#     #         neighbor_col.append(neighbor_col_value)
#     #         neighbor_col_value+=1
#     #     neighbor_col_value=-around_grid_num

#     # for i in range(around_grid_num):
#     #     for j in range(around_grid_num+1):
#     #         neighbor_col.append(neighbor_col_value)
#     #         neighbor_col_value+=1
#     #     neighbor_col_value=-around_grid_num

#     return neighbor_row,neighbor_col
    

# def make_pairlist_g(pair_list_g,select_gx,select_gy,G_MAP,select_index,NEIGHBOR_LEN,NEIGHBOR_COL,NEIGHBOR_ROW):
#     particle_counter=0
#     for k in range(NEIGHBOR_LEN):
#         search_gx=select_gx+NEIGHBOR_COL[k]
#         search_gy=select_gy+NEIGHBOR_ROW[k]
#         if(search_gx>=N_GX):
#             search_gx-=N_GX
#         elif(search_gx<0):
#             search_gx+=N_GX
#         if (search_gy>=N_GY):
#             search_gy-=N_GY
#         elif(search_gy<0):
#             search_gy+=N_GY
        
#         search_index=G_MAP[search_gy][search_gx]
#         pair_list_g[select_index][particle_counter]=search_index

#         particle_counter+=1

#     pair_list_g[select_index][-1]=particle_counter





    


# def gmap_create(l_gx,l_gy,pair_length,NP,N_GX,N_GY,rx,ry,neighbor_row,neighbor_col):
#     G_MAP=np.full((N_GY,N_GX),-1)
#     pair_list_g=np.full(pair_length,-1)
#     for i in range(NP):
#         gx_map=int(rx[i]/l_gx)
#         gy_map=int(ry[i]/l_gy)
#         G_MAP[gy_map][gx_map]=i
#     for i in range(N_GY):
#         for j in range(N_GX):
#             select_index=G_MAP[i][j]
#             make_pairlist_g()

#     return pair_list_g
