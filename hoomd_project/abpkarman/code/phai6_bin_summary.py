import numpy as np 
import matplotlib.pyplot as plt
import math
import os 
import numpy.polynomial.polynomial as P
path="./"
files=os.listdir(path)
files_dir=[f for f in files if os.path.isdir(os.path.join(path,f))]
print(files_dir)

# fiels_dir=['0.8_10.0_40.0_1.5_0.1', '0.8_10.0_40.0_8.0_0.5', '0.8_3.0_40.0_0.5_0.1', '0.8_5.0_40.0_10.0_0.5',
#   '0.8_5.0_40.0_8.0_0.1', '0.8_5.0_40.0_2.0_0.1', '0.8_3.0_40.0_6.0_0.5', '0.8_5.0_40.0_4.0_0.1', 
#   '0.8_10.0_40.0_4.0_0.5', '0.8_5.0_40.0_6.0_0.5', '0.8_3.0_40.0_2.5_0.1', '0.8_5.0_40.0_0.0_0.1', '0.8_5.0_40.0_5.0_0.5', '0.8_3.0_40.0_2.0_0.1', 
#  '0.8_10.0_40.0_5.0_0.5', '0.8_3.0_40.0_3.0_0.1', '0.8_3.0_40.0_0.0_0.5', '0.8_3.0_40.0_6.0_0.1', 
#  '0.8_5.0_40.0_4.0_0.5', '0.8_10.0_40.0_10.0_0.5', '0.8_5.0_40.0_1.5_0.1', '0.8_10.0_40.0_2.0_0.5', 
#  '0.8_10.0_40.0_0.0_0.5', '0.8_3.0_40.0_10.0_0.5', '0.8_10.0_40.0_0.0_0.1', '0.8_5.0_40.0_2.5_0.1',
#    '0.8_10.0_40.0_1.0_0.1', '0.8_3.0_40.0_0.0_0.1', '0.8_3.0_40.0_8.0_0.5', '0.8_5.0_40.0_7.0_0.5', 
#    '0.8_10.0_40.0_0.5_0.1', '0.8_3.0_40.0_8.0_0.1', '0.8_5.0_40.0_0.5_0.1', '0.8_10.0_40.0_7.0_0.5', 
#    '0.8_10.0_40.0_6.0_0.1', '0.8_3.0_40.0_5.0_0.5', '0.8_5.0_40.0_1.0_0.1', '0.8_10.0_40.0_2.0_0.1', 
#    '0.8_3.0_40.0_1.0_0.1', '0.8_3.0_40.0_1.5_0.1', '0.8_5.0_40.0_3.0_0.1', '0.8_3.0_40.0_4.0_0.5',
#      '0.8_10.0_40.0_4.0_0.1', '0.8_10.0_40.0_3.0_0.1', '0.8_3.0_40.0_4.0_0.1', '0.8_3.0_40.0_2.0_0.5', 
#      '0.8_10.0_40.0_8.0_0.1', '0.8_5.0_40.0_0.0_0.5', '0.8_5.0_40.0_6.0_0.1', '0.8_3.0_40.0_7.0_0.5', 
#      '0.8_5.0_40.0_2.0_0.5', '0.8_10.0_40.0_2.5_0.1', '0.8_10.0_40.0_6.0_0.5', '0.8_5.0_40.0_8.0_0.5']



# files_dir=[f for f in  files_dir if "0.5" not in f]
# print(files_dir)

rho=0.8 
ave_flow_list=[5.0]
static_dia=50.0
acive_force_list=[0.0 , 1.0  , 5.0,10.0] 

DR=10.0


files_dir=[]
for ave_flow in ave_flow_list:
    for active_force in acive_force_list:

        temp_dir="{0}_{1}_{2}_{3}_{4}".format(rho,ave_flow,static_dia,active_force,DR)
        files_dir.append(temp_dir)
print(files_dir)
print(len(files_dir))
lx=static_dia*14
ly=static_dia*8.0

static_rx=lx/4
static_ry=ly/2
# 検査体積 binn
l_grid_bin=1.5
n_grid_bin=math.ceil(lx/l_grid_bin)
Vc=l_grid_bin*ly
bin_x=[i*l_grid_bin for i in range(n_grid_bin)]
print(len(bin_x))

plt.figure(figsize=(7.0 ,4.0),dpi=300)
for dir in files_dir:
    try:
        data=np.load(dir+"/phi6_bin2.npy")
    except:
        print("skip".format(dir))
        continue


    print(dir)
    # print(len(data))
    # data=data[9:]
    print(data.shape)
    mean_data=np.mean(data,axis=0)
    print(mean_data.shape)
    coefs=np.polyfit(bin_x,mean_data,10)
    
    # print(coefs)
    fitted_curve=np.poly1d(coefs)(bin_x)
    # print(fitted_curve)
    v0=dir.split("_")[3]
    U=dir.split("_")[1]
    plt.ylim(0.4 ,0.95)
    plt.scatter(bin_x[0::1],mean_data[0::1],label=r'$v_0={0}$'.format(v0))
    plt.axvspan(static_rx-static_dia/2,static_rx+static_dia/2,color="blue",alpha=0.05)
    # plt.plot(bin_x,fitted_curve,label=dir)
    # output_path=dir+"/phi6_summary.png"
    plt.title(r'$U={0},DR={1}$'.format(U,DR))
    plt.xlabel(r"$x$",fontsize=15)
    plt.ylabel(r'$\|\phi_6\|$',fontsize=15)

output_path="./phi6_summary2_U={0}_DR={1}.png".format(ave_flow_list[0],DR)

plt.legend(loc="upper left",bbox_to_anchor=(1,1.0))
plt.savefig(output_path,bbox_inches="tight")

    # plt.cla()
    # plt.clf()
    