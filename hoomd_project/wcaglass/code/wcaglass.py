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
num_out=int(8000)
# 出力をしているときのステップ数の間隔
pos_out_steps_period=int(pos_out_steps/num_out)
##############################################################
# image出力枚数
image_out_num=150
# 総出力回数に対してイメージを出力させる間隔
image_out_period=int(num_out/image_out_num)




# pos_hout=int(nsteps/500)
# thermo_hout=int(nsteps/2.0)

data_file="../data/nu0.8.npz"

npz=np.load(data_file)

rx=npz["rx"]
ry=npz["ry"]
rx=rx-lx/2.0
ry=ry-ly/2.0
# sigma=npz["sigma"]
# typeid=npz["typeid"]
print(N)



temp_radius = np.sqrt(lx*ly*nu/((N0+N1*AL**2)*np.pi)) #radius

def binary_a(radius,N,ly,N1,AL):
    sigma = [0 for i in range(N)]
    N1_num=0
    for i in range(0,N-1,2):
        K = int(i/ly)+1
        if K%2==1:
            sigma[i] = radius
            if N1_num<=N1:
                sigma[i+1]=AL*radius
                N1_num+=1
            else:
                sigma[i+1]=radius
        
        else:
            if N1_num<=N1:
                sigma[i]=AL*radius
                N1_num+=1
            else:
                sigma[i]=radius
            
            sigma[i+1]=radius
    return sigma
         


sigma=binary_a(temp_radius,N,ly,N1,AL)
print(len(sigma))
sigma=[sigma[i]*2 for i in range(N)]
print(len(sigma))



set_sigma=list(set(sigma))
small_sigma=set_sigma[0]
large_sigma=set_sigma[1]

typeid=[0 for i in range(N)]
for i in range(N):
    if sigma[i]==small_sigma:
        typeid[i]=0
    else:
        typeid[i]=1


sigma_ss=small_sigma
sigma_ll=large_sigma
sigma_sl=(small_sigma+large_sigma)/2



ver=str(nu)+"_"+str(fixed_percent)+"_"+str(kbT)

main_dir="./"+ver
if not os.path.exists(main_dir): os.makedirs(main_dir)
print(main_dir)





fixed_id=[]

fixed_N=int(N*fixed_percent)
move_N=N-fixed_N


while len(fixed_id)<fixed_N:
    id=random.randint(0,N-1)
    if id not in fixed_id:
        fixed_id.append(id)

move_id=[num for num in range(N) if num not in fixed_id]

print(len(move_id))

position=[]
for i in range(len(rx)):
    position.append((rx[i],ry[i],0.0))



snapshot = gsd.hoomd.Frame()
snapshot.particles.N = N
snapshot.particles.position = position[0:N]



snapshot.particles.typeid = typeid
snapshot.particles.types = ['small','large']
snapshot.particles.diameter=sigma
snapshot.particles.mass = np.ones((N))
snapshot.configuration.box = [lx, ly, 0, 0, 0, 0]


sim = hoomd.Simulation(device=hoomd.device.GPU(), seed=12)

sim.create_state_from_snapshot(snapshot)
# Integration information


cell = hoomd.md.nlist.Cell(buffer=0.4)
lj = hoomd.md.pair.LJ(nlist=cell, mode="shift")
lj.params[("small", "small")] = dict(epsilon=epsilon, sigma=sigma_ss)
lj.r_cut[("small", "small")] = 2**(1/6)*sigma_ss


lj.params[("small", "large")] = dict(epsilon=epsilon, sigma=sigma_sl)
lj.r_cut[("small", "large")] = 2**(1/6)*sigma_sl

lj.params[("large", "large")] = dict(epsilon=epsilon, sigma=sigma_ll)
lj.r_cut[("large", "large")] = 2**(1/6)*sigma_ll


nve = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.Tags(move_id))
integrator = hoomd.md.Integrator(
 dt=dt,
 methods=[nve],
 forces=[lj],
)

# velocity_operation=hoomd.update.CustomUpdater(action=RelativeFlow(ave_flow,dt),trigger=1)
# sim.operations+=velocity_operation
# sim.operations += rotational_diffusion_updater
sim.operations.integrator = integrator

# カスタムクラス
class PrintTimestep(hoomd.custom.Action):
    def act (self,timestep):
        print(timestep)
custom_action = PrintTimestep()
custom_op = hoomd.write.CustomWriter(action=custom_action,
                                 trigger=hoomd.trigger.Periodic(10000))
sim.operations.writers.append(custom_op)





# logger定義
pos_logger = hoomd.logging.Logger()
# pos_logger.add(sim, quantities=['timestep', 'walltime'])
# pos_logger.add(thermodynamic_properties,quantities=["kinetic_temperature","kinetic_energy","potential_energy","volume"])
gsd_writer_pos = hoomd.write.GSD(filename=main_dir+"/log_pos.gsd",
                            #  trigger=hoomd.trigger.Periodic(pos_hout),
                            trigger=hoomd.trigger.And([
                                    hoomd.trigger.After(after_steps),
                                    hoomd.trigger.Periodic(pos_out_steps_period)]),
                             mode='xb',
                             filter=hoomd.filter.All())
gsd_writer_pos.log = pos_logger

sim.operations.writers.append(gsd_writer_pos)

sim.state.thermalize_particle_momenta(filter=hoomd.filter.Tags(move_id), kT=kbT)

sim.run(0)

print("-----run--------")
second_time=time.time()



variables = {
    'real_time': real_time,
    'after_time': after_time,
    'pos_out_time': pos_out_time,
    'dt': dt,
    'nsteps': nsteps,
    'after_steps': after_steps,
    'pos_out_steps': pos_out_steps,
    'pos_out_steps_period': pos_out_steps_period,
    'num_out': num_out,
    'image_out_num': image_out_num,
    'image_out_period': image_out_period,
    'm': m,
    'epsilon': epsilon,
    'kbT': kbT,
    'nu': nu,
    'fixed_percent': fixed_percent,
    'sigma_ss': sigma_ss,
    'sigma_ll': sigma_ll,
    'sigma_sl': sigma_sl
}

# ファイルにほぞん.txt
with open(main_dir+'/data.txt', 'w') as f:
    for var_name, var_value in variables.items():
        f.write(f'{var_name} {var_value}\n')


# ファイルに保存.pickle
with open(main_dir+'/data.pickle', 'wb') as handle:
    pickle.dump(variables, handle, protocol=pickle.HIGHEST_PROTOCOL)

# # ファイルから読み込み
# with open('data.pickle', 'rb') as handle:
#     loaded_dict = pickle.load(handle)

sim.run(nsteps)

print(time.time()-first_time)

print(time.time()-second_time)

traj = gsd.hoomd.open(main_dir+"/log_pos.gsd",mode='r')


output_dir=main_dir+"/figure_2d"
if not os.path.exists(output_dir): os.makedirs(output_dir)
figsize=(10,10)
plt.figure(figsize=figsize,dpi=150)

print(len(traj))

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
    
