import os 
from PIL import Image
import sys
# ver=1
rho=float(sys.argv[1])
V0=float(sys.argv[2])
# static_dia=float(sys.argv[3])
# red_v=float(sys.argv[4])
# rotate_dif=float(sys.argv[5])
# # N=2500

# ly=static_dia*5
# lx=static_dia*18
# red_v=format(red_v,'.6f')
# rotate_dif=format(rotate_dif,'.6f')
# rho=format(rho,'.2f')
# va=format(va,'.2f')
# ver="rho_{:.2f}_v0_{:.1f}".format(temp_rho,V0)
ver="test"
main_dir="./{0}/".format(ver)
figure_dir=main_dir+"figure/"

# figure_dir=main_dir+"figure/"

# main_dir="./"+str(rho)+"_"+str(ave_flow)+"_"+str(static_dia)+"_"+str(red_v)+"_"+str(rotate_dif)+"/"

############アニメーション################    
images=[]
# image_num=sum(os.path.isfile(os.path.join(pic_output name)) for name in os.listdir(pic_output))
image_num=sum(os.path.isfile(os.path.join(figure_dir,name))for name in os.listdir(figure_dir))
print(image_num)
for i in range(1,image_num,2):
    file_name=figure_dir+"/figure"+str(i)+".png"
    im=Image.open(file_name)
    images.append(im)
# if not os.path.exists(output_dir): os.makedirs(output_dir)
images[0].save(main_dir+"out_ela2.gif",save_all=True,append_images=images[1:],loop=0,duration=10)
    