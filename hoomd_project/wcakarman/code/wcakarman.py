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
def initialize_ra(Nx,Ny,dx,dy,lx,ly,remove_dx_num):
    rx=[]
    ry=[]
#     rx.append(lx/4.0-lx/2)
#     ry.append(ly/2.0-ly/2)
    
    for i in range(Nx):
        for j in range(Ny):
            if (j>=int(Ny/2)-remove_dx_num and j<=int(Ny/2)+remove_dx_num):
                if(i>=int(Nx/4)-remove_dx_num and i<=int(Nx/4) +remove_dx_num):
                    continue


            # rx.append(i*dx+dx/2)
            # ry.append(j*dy+dy/2)
            
            rx.append(i*dx-lx/2+dx/2)
            ry.append(j*dy-ly/2+dy/2)



    return rx,ry

class RelativeFlow(hoomd.custom.Action):
    def __init__(self, ave_flow,h):
        self.ave_flow=ave_flow
        self.h=h
        
    def act (self,timestep):
        snap=self._state.get_snapshot()
        if snap.communicator.rank==0:
            snap.particles.position[:,0]+=self.ave_flow*self.h
        self._state.set_snapshot(snap)
# class ChangeVelocity(hoomd.custom.Action):
#     def __init__(self, velocity,lx):
#         self.velocity=velocity
#         self.lx=lx
        
#     def act (self,timestep):
#         snap=self._state.get_snapshot()
#         if snap.communicator.rank==0:
#             position=snap.particles.position
#             position_x=position.T[0]
#             # change_index=np.where(position_x>lx-1)
#             # snap.particles.velocity[change_index,0]=self.velocity
#             snap.particles.velocity[:,0]+=self.velocity
#         self._state.set_snapshot(snap)
rho=float(sys.argv[1])
ave_flow=float(sys.argv[2])
static_dia=float(sys.argv[3])
reduced_speed=float(sys.argv[4])
rotational_diffusion=float(sys.argv[5])

pattern=str(rho)+"_"+str(ave_flow)+"_"+str(static_dia)+"_"+str(reduced_speed)+"_"+str(rotational_diffusion)

init_pattern=str(rho)+"_"+str(0.0)+"_"+str(static_dia)+"_"+str(0.0)+"_"+str(5.0)

version="ver1"

main_dir="./"+pattern




input_traj = gsd.hoomd.open('../../data/pattern4/{0}/log_pos.gsd'.format(init_pattern))

input_pos=input_traj[-1].particles.position.T
rx=input_pos[0]
ry=input_pos[1]


dx=np.sqrt(1.0/rho)
dy=dx
lx=input_traj[0].configuration.box[0]
ly=input_traj[0].configuration.box[1]
Nx=int(lx/dx)
Ny=int(ly/dy)


figsize=(10.0,10.0*ly/lx)
remove_dx_num=math.ceil(static_dia*0.5/dx)


# temp=1.0
m=1.0
epsilon=1.0
kbT = 1.0



# 刻み幅小さすぎの可能性もあるから大きめにしてみてもいいかも

# h=5e-4
# dt = 1e-6
real_time=250
dt = 1e-4
nsteps=int(real_time/dt)
pos_hout=int(nsteps/500)

# t_end=nsteps*h
# h_rev=1/h 
# h2=0.5*h*h



# rx,ry=initialize_ra(Nx,Ny,dx,dy,lx,ly,remove_dx_num)
# rx=np.loadtxt("../data_file/square/rx_{0}_{1}.txt".format(rho,static_dia))
# ry=np.loadtxt("../data_file/square/ry_{0}_{1}.txt".format(rho,static_dia))

N=len(rx)
print(N)
position=[]
for i in range(len(rx)):
    position.append((rx[i],ry[i],0.0))



snapshot = gsd.hoomd.Frame()
snapshot.particles.N = N
snapshot.particles.position = position[0:N]


# (1,0,0,0,)
# snapshot.particles.orientation = orientation
snapshot.particles.typeid = (N)*[0]
snapshot.particles.types = ['Move']
snapshot.particles.diameter=(N)*[1.0]
snapshot.particles.mass = np.ones((N))
snapshot.configuration.box = [lx, ly, 0, 0, 0, 0]


# with gsd.hoomd.open(name='_init'+str(static_dia)+'.gsd', mode='xb') as f:
#     f.append(snapshot)





print(main_dir)
output_dir=main_dir+"/figure_2d"
if not os.path.exists(output_dir): os.makedirs(output_dir)

sim = hoomd.Simulation(device=hoomd.device.GPU(), seed=12)

sim.create_state_from_snapshot(snapshot)
# Integration information


sigma_aa=1.0
sigma_ab=(static_dia/2+0.5)
epsilon_aa=1.0
epsilon_ab=1.0

cell = hoomd.md.nlist.Cell(buffer=0.4)
lj = hoomd.md.pair.LJ(nlist=cell, mode="shift")
lj.params[("Move", "Move")] = dict(epsilon=epsilon_aa, sigma=sigma_aa)
lj.r_cut[("Move", "Move")] = 2**(1/6)*sigma_aa


walls=[hoomd.wall.Cylinder(origin=(lx/4-lx/2,ly/2-ly/2,0),radius=static_dia/2,inside=False,axis=(0,0,1))]
ljw = hoomd.md.external.wall.LJ(walls=walls)
# ljw.params['Move'] = {"epsilon": epsilon_ab, "sigma": sigma_ab, "r_cut": sigma_ab*2**(1/6)}

ljw.params['Move'] = {"epsilon": 1.0, "sigma": sigma_aa/2.0, "r_cut": sigma_ab*2**(1/6)-static_dia/2}



nve = hoomd.md.methods.ConstantVolume(filter=hoomd.filter.All())
integrator = hoomd.md.Integrator(
 dt=dt,
 methods=[nve],
 forces=[lj, ljw],
)

velocity_operation=hoomd.update.CustomUpdater(action=RelativeFlow(ave_flow,dt),trigger=1)
sim.operations+=velocity_operation
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
# pos_logger.add(thermodynamic_properties,quantities=["kinetic_temperature","kinetic_energy","potential_energy","volume"])
gsd_writer_pos = hoomd.write.GSD(filename=main_dir+"/log_pos.gsd",
                             trigger=hoomd.trigger.Periodic(pos_hout),
                             mode='xb',
                             filter=hoomd.filter.All())
gsd_writer_pos.log = pos_logger

sim.operations.writers.append(gsd_writer_pos)


sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kbT)
# 初期速度変化
sim.run(0)




variables = {
    'real_time': real_time,
    'dt': dt,
    'nsteps': nsteps,
    'N': N,
    'm': m,
    'epsilon_aa': epsilon_aa,
    'epsilon_ab': epsilon_ab,
    
    'kbT': kbT,
    'rho': rho,
    'ave_flow': ave_flow,
    'static_dia': static_dia,
    'sigma_aa': sigma_aa,
    'sigma_ab': sigma_ab,
}



# ファイルにほぞん.txt
with open(main_dir+'/data.txt', 'w') as f:
    for var_name, var_value in variables.items():
        f.write(f'{var_name} {var_value}\n')


# ファイルに保存.pickle
with open(main_dir+'/data.pickle', 'wb') as handle:
    pickle.dump(variables, handle, protocol=pickle.HIGHEST_PROTOCOL)






print("-----run--------")
second_time=time.time()

sim.run(nsteps)

print(time.time()-first_time)

print(time.time()-second_time)
# os.chdir("../")
# traj = gsd.hoomd.open('./'+ver+'/log_pos_'+ver+'.gsd', 'rb')


# # traj = gsd.hoomd.open('log_force2d_'+ver+'.gsd', 'rb')
# plt.figure(figsize=figsize)

# print(len(traj))

# for t in range(len(traj)-1,0,-2):
#     print(t)
#     bx=plt.axes()
#     plt.axis([-lx/2,lx/2,-ly/2,ly/2])
#     position=traj[t].particles.position   
#     c=pat.Circle(xy=(lx/4-lx/2,ly/2-ly/2),radius=static_dia/2,fc="b")
#     bx.add_patch(c) 
#     for i in range(N):
#         c=pat.Circle(xy=(position[i][0],position[i][1]),radius=0.5,fc="r")
#         bx.add_patch(c)

#     plt.title("step"+str(t))
#     plt.savefig(output_dir+"/figure{0}.png".format(t))
#     plt.cla()

#     ############アニメーション################    
# images=[]
# # image_num=sum(os.path.isfile(os.path.join(pic_output name)) for name in os.listdir(pic_output))
# image_num=sum(os.path.isfile(os.path.join(output_dir,name))for name in os.listdir(output_dir))
# print(image_num)
# for i in range(0,image_num):
#     file_name=output_dir+"/figure"+str(i)+".png"
#     im=Image.open(file_name)
#     images.append(im)

# gif_output_dir=main_dir+"/abpgif2"

# if not os.path.exists(gif_output_dir): os.makedirs(gif_output_dir)
# images[0].save(gif_output_dir+"/out_ela2.gif",save_all=True,append_images=images[1:],loop=0,duration=10)
    
    