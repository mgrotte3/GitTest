# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 22:33:05 2016

@author: mario90
"""
import numpy as np
# Set up matplotlib and use a nicer set of plot parameters
#%config InlineBackend.rc = {}
#import matplotlib
#matplotlib.rc_file("../../templates/matplotlibrc")
import matplotlib.pyplot as plt
import scipy.stats
#%matplotlib inline
from astropy.io import fits

def image_statistics(data):
    print('Min:', np.min(data))
    print('Max:', np.max(data))
    print('Mean:', np.mean(data))
    print('Stdev:', np.std(data))
    print(type(data.flat))
    
image_list = [ 'C:\\Users\\mario90\\Desktop\\AllSky ISS\\AllSkyImage00018289'+n+'.FIT' \
              for n in ['3','4','5','6'] ]

# The long way
image_concat = []
for image in image_list:
    image_concat.append(fits.getdata(image))

# The short way
#image_concat = [ fits.getdata(image) for image in IMAGE_LIST ]

# The long way
final_image = np.zeros(shape=image_concat[0].shape)

for image in image_concat:
    final_image += image

# The short way
#final_image = np.sum(image_concat, axis=0)
imag1 = image_concat[-1]
cropped_image_data1=imag1[381:418,206:243]
cropped_image_data2=final_image[381:418,206:243]
plt.figure(1)
image_hist = plt.hist(final_image.flat, 1000)
plt.figure(2)
plt.imshow(final_image, cmap='gray', vmin=3.e4, vmax=8.e4)
plt.colorbar()
plt.figure(3)
plt.imshow(cropped_image_data1, cmap='gray')
plt.colorbar()
image_statistics(cropped_image_data1)
plt.figure(4)
plt.imshow(cropped_image_data2, cmap='gray')
plt.colorbar()
image_statistics(cropped_image_data2)
snr1=scipy.stats.signaltonoise(cropped_image_data1, axis=None)
snr2=scipy.stats.signaltonoise(cropped_image_data2, axis=None)
print('SNR of Original Image:', snr1)
print('SNR of Stacked Image:', snr2)