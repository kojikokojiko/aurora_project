
import matplotlib.pyplot as plt
import matplotlib.patches as pat

import math

rx = []
ry = []
def initialize_r(rx, ry, static_dia, lx, ly, temp_rho):
    dx = math.sqrt(1.0 / temp_rho)
    dy = dx
    nx = int(lx / dx)
    ny = int(ly / dy)

    remove_dx_num = math.ceil(static_dia * 0.5 / dx)
    for i in range(nx):
        for j in range(ny):
            x = i * dx + dx / 2 - 0.5 * lx
            y = j * dy + dy / 2 - 0.5 * ly
            if x < 0:
                r2 = x * x + y * y
                if r2 < (0.5*static_dia + 0.5)*(0.5*static_dia + 0.5) :
                    continue
            else:
                if y < static_dia / 2.0 and y > -static_dia / 2.0:
                    if x < 0.5:
                        continue
            rx.append(x)
            ry.append(y)

# Usage:

static_dia=100
L=200
lx=L
ly=L
temp_rho=0.8
initialize_r(rx, ry, static_dia, lx, ly, temp_rho)


N=len(rx)
print(rx)
print(N)
plt.figure(figsize=(6,6))
bx=plt.axes()
plt.axis([-lx/2,lx/2,-ly/2,ly/2])
semicircle = pat.Wedge(center=(0, 0), r=static_dia/2, theta1=90, theta2=270, facecolor=[0,0,1,0.5])
bx.add_patch(semicircle)

for i in range(N):
    c=pat.Circle(xy=(rx[i],ry[i]),radius=0.5,fc="r")
    bx.add_patch(c)


# plt.title("step"+str(t))
plt.savefig("./figure{0}.png".format(N))
plt.cla()
