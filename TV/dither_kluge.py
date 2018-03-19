from pylab import *
from PIL import Image
from scipy import *
from FFT import *

inputfilename='/projects/souris/jgleave/TissueVision_2012/23Oct22012/DCXz22.tif'
outputfilename='/projects/souris/bjnieman/JG_files/stains/DCXz22_23Oct2012.tif'

imgfh = Image.open(inputfilename)
img = reshape( array(imgfh.getdata(),float) ,(imgfh.size[0],imgfh.size[1]))

imgorig = img.copy()
figure(1); imshow(imgorig)

for j in range(16):
    eopshift=-1.5+0.25*j
    img = imgorig.copy()
    phaseramp = exp(-arange(img.shape[1],dtype=float)*1.j*2*pi*eopshift/img.shape[1])
    img[0::2,:]=abs( fftshift(ifft(fftshift( fftshift(fft(fftshift(img[0::2,:],axes=(1,)),axis=1),axes=(1,))\
                                             *phaseramp[newaxis,:],axes=(1,)),axis=1),axes=(1,)) )
    figure(j+1); imshow(img)

#shift is clearly not uniform over the tile!!
#need to include positional dependence for the shift

#emprical play
pos=array([14,29,75,114,150,213,228,234],float)     #xpos
shift=array([0,0.5,0,-0.9,-0.5,1.3,0.5,-0.5],float) #req'd shift
spfit=UnivariateSpline(pos,shift,s=0)
fig1=figure(1); plot(pos,shift,'o'); plot(arange(img.shape[1]),spfit(arange(img.shape[1])),'-')

#try variable shift as prescribed by emprical play above
modpos=arange(img.shape[1])+spfit(arange(img.shape[1]))

img=imgorig.copy()
for j in range(0,img.shape[0],2):
    f=interp1d(modpos,img[j,:],kind='cubic',bounds_error=False,fill_value=0.0)
    img[j,:]=f(arange(img.shape[1]))

modimgfh=Image.fromarray((255.0*img/max(img.flat)).astype(int8),mode="L")
modimgfh.save(outputfilename)

