import matplotlib.pyplot as plt 
import matplotlib.image as mpimg 
import numpy as np
from PIL import Image

probe = np.array(Image.open('probe1.png'))
probe0 = Image.open('probe1.png')
plt.imshow(probe)
pos1=plt.ginput(2)
print(pos1)
distance1 = abs(pos1[0][1]-pos1[1][1])
ytip = max(pos1[0][1],pos1[1][1])
ysurf = min(pos1[0][1],pos1[1][1])
ephys = np.array(Image.open('scale.png'))
ephys0 = Image.open('scale.png')
plt.imshow(ephys)
pos2=plt.ginput(2)
print(pos2)
distance2 = abs(pos2[0][1]-pos2[1][1])
crop_probe = probe0.crop((0,50,150,1012))
sizex = 150*distance2/distance1
sizey = 962*distance2/distance1
probeResize = crop_probe.resize((int(sizex),int(sizey)))
result = Image.new('RGB',(1800,3000))
result.paste(ephys0,(0,int((ysurf-50)*distance2/distance1)))
result.paste(probeResize,(1314,0))
plt.imshow(result)
plt.show()
