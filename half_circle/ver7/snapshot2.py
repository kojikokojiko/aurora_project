import matplotlib.pyplot as plt
import matplotlib.patches as pat
import numpy as np
import os
import sys

temp_rho=float(sys.argv[1])
V0=float(sys.argv[2])
ver="rho_{:.2f}_v0_{:.1f}".format(temp_rho,V0)
print(ver)
main_dir="./{0}/".format(ver)
figure_dir=main_dir+"figure/"
if not os.path.exists(figure_dir):
    os.makedirs(figure_dir)
rx_list=np.loadtxt(main_dir+"rx.dat")
ry_list=np.loadtxt(main_dir+"ry.dat")
print(rx_list.shape)

N=rx_list.shape[1]
static_dia=50
lx=static_dia*2
ly=static_dia*2 


plt.figure(figsize=(6,6))
bx=plt.axes()

for t in range(len(rx_list)):
    print(t)
    rx=rx_list[t]
    ry=ry_list[t]
        

    plt.axis([-lx/2,lx/2,-ly/2,ly/2])
    semicircle = pat.Wedge(center=(0, 0), r=static_dia/2, theta1=90, theta2=270, facecolor=[0,0,1,0.5])
    bx.add_patch(semicircle)
    
    r_c_semicircle = pat.Wedge(center=(0, 0), r=(static_dia+1.0)/2*pow(2,1.0/6.0), theta1=90, theta2=270, facecolor=[0,1,0,0.4])
    bx.add_patch(r_c_semicircle)
    
    
    for i in range(N):
        c=pat.Circle(xy=(rx[i],ry[i]),radius=0.5,fc="r")
        bx.add_patch(c)


    plt.title("step"+str(t))
    plt.savefig(figure_dir+"figure{0}_test.png".format(t))
    plt.cla()
