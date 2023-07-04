import numpy as np
import matplotlib.patches as pat
import matplotlib.pyplot as plt

def initialize_ra(Nx,Ny,dx,dy,lx,ly):
    rx=[]
    ry=[]
#     rx.append(lx/4.0-lx/2)
#     ry.append(ly/2.0-ly/2)
    
    for i in range(Nx):
        for j in range(Ny):    
            rx.append(i*dx-lx/2+dx/2)
            ry.append(j*dy-ly/2+dy/2)
    return rx,ry

rho=0.8
static_dia=50.0
dx=np.sqrt(1.0/rho)
dy=dx
lx=static_dia*14
ly=static_dia*8
Nx=int(lx/dx)
Ny=int(ly/dy)
sigma=1.0
static_rx=lx/4.0-lx/2.0
static_ry=ly/2.0-ly/2.0



temp_rx,temp_ry=initialize_ra(Nx,Ny,dx,dy,lx,ly)

temp_NP=len(temp_rx)

fine_index=[]
rc=(static_dia+sigma)/2.0
rc2=rc*rc
for i in range(temp_NP):
    rxij=static_rx-temp_rx[i]
    ryij=static_ry-temp_ry[i]

    r2=rxij*rxij+ryij*ryij

    if(r2>rc2):
        fine_index.append(i)


rx= [temp_rx[i] for i in fine_index]
ry= [temp_ry[i] for i in fine_index]


NP=len(rx)

print(NP)




# 描画設定##############################################
# https://qiita.com/stanaka2/items/c40841f858d7083aad4e
size=1.0
ppi=72 # points per inche 

fig = plt.figure(figsize=(14,8), dpi=200.0)
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

# plt.figure(figsize=(12,8))
# bx=plt.axes()
# plt.axis([-lx/2,lx/2,-ly/2,ly/2])

 
# for i in range(NP):
#     c=pat.Circle(xy=(rx[i],ry[i]),radius=0.5,fc="r")
#     bx.add_patch(c)
im=ax.scatter(rx,ry,s=size**2,c="r", linewidths=0)
c=pat.Circle(xy=(lx/4-lx/2,ly/2-ly/2),radius=static_dia/2,fc="b")
ax.add_patch(c)
plt.title("step")
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

ax.set_xticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 5))
ax.set_yticks(np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 5))

plt.savefig("./figure{0}_{1}.png".format(rho,static_dia))
plt.cla()


rx=np.array(rx)
ry=np.array(ry)
np.savetxt("./rx_{0}_{1}.txt".format(rho,static_dia),rx)
np.savetxt("./ry_{0}_{1}.txt".format(rho,static_dia),ry)

# np.savetxt("./NNnum{0}.txt".format(t),NNnum)

