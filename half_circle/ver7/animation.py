



from PIL import Image


# added
import os
import sys

rho=float(sys.argv[1])
V0=float(sys.argv[2])
Image.MAX_IMAGE_PIXELS = 1000000000

ver=5
nsteps=2*10**4
main_dir="./rho_{:.2f}_v0_{:.1f}/".format(rho,V0)
output_dir=main_dir+"gif/"
pic_output=main_dir+"figure/"


# 画像ファイルのリストを取得
image_files = os.listdir(pic_output)
# image_files.sort()  # ファイル名のソート
len_files=len(image_files)
print(image_files)
###########アニメーション################    
images=[]

for i in range(len_files):
    file_name="figure{0}_test.png".format(i)

    file_path = os.path.join(pic_output, file_name)
    im=Image.open(file_path)
    images.append(im)
if not os.path.exists(output_dir): os.makedirs(output_dir)
images[0].save(output_dir+"/out_ela.gif",save_all=True,append_images=images[1:],loop=0,duration=30)
    
        
    
