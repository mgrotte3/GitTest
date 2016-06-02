# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 23:07:21 2016

@author: mario90
"""
import pywt
import numpy as np
# Set up matplotlib and use a nicer set of plot parameters
#%config InlineBackend.rc = {}
#import matplotlib
#matplotlib.rc_file("../../templates/matplotlibrc")
import matplotlib.pyplot as plt
import scipy.stats
#%matplotlib inline
from astropy.io import fits

def wavelet_fusion(Im1, Im2):
    m1 = np.mean(Im1)
    m2 = np.mean(Im2)
    s1 = np.std(Im1)
    s2 = np.std(Im2)
    g = s2/s1
    offset = m2-g*m1
    ImH = Im1*g+offset
    # Wavelet Transform step
    coeffs1 = pywt.dwt2(ImH, 'db1')
    coeffs2 = pywt.dwt2(Im2, 'db1')
    # SELECTION OF COEFFICIENT USING Linear combination BETWEEN both approximations.
    A = 2*(coeffs1[0]+coeffs2[0])/2
    # INVERSE WAVELET TRANSFORM TO GENERATE THE FUSED IMAGE.
    Xsyn = pywt.idwt2((A,coeffs1[1]), 'db1')
    return Xsyn
    
if __name__ == "__main__":
    image1 = 'C:\\Users\\mario90\\Desktop\\AllSky ISS\\AllSkyImage000182894.FIT'
    image2 = 'C:\\Users\\mario90\\Desktop\\AllSky ISS\\AllSkyImage000182895.FIT'
    image3 = 'C:\\Users\\mario90\\Desktop\\AllSky ISS\\AllSkyImage000182896.FIT'
    imag1 = fits.getdata(image1)
    imag2 = fits.getdata(image2)
    imag3 = fits.getdata(image3)
    Xsyn_inter = wavelet_fusion(imag1,imag2)
    Xsyn_final = wavelet_fusion(Xsyn_inter, imag3)
    cropped_image_data1=imag3[381:418,206:243]
    cropped_image_data2=Xsyn_final[381:418,206:243]
    plt.figure(1)
    plt.imshow(Xsyn_final, cmap='gray')
    plt.colorbar()
    plt.figure(2)
    plt.imshow(cropped_image_data1, cmap='gray')
    plt.colorbar()
    plt.figure(3)
    plt.imshow(cropped_image_data2, cmap='gray')
    plt.colorbar()
    snr1 = scipy.stats.signaltonoise(cropped_image_data1, axis=None)
    snr2 = scipy.stats.signaltonoise(cropped_image_data2, axis=None)
    print('SNR of Original Image:', snr1)
    print('SNR of Stacked Image:', snr2)
#ima1=np.max(imag2)
#imi1=np.min(imag2)
#mse1=np.std(imag2)
#snr1=20*np.log10((ima1-imi1)/mse1)
##fprintf(1, 'SNR of Original Image: %f\n',snr1)
##fprintf(1, 'Original: Number of Pixels: %f\n',numel(im2bw(X2(:))))
##fprintf(1, 'Original: Number of True Pixels: %f\n',sum(im2bw(X2(:))))
##fprintf(1, 'Original: Maximum Pixel Value relative to Orginial: %f\n',max(X2(:))) % Max Pixel value
#ima2=np.max(Xsyn_final)
#imi2=np.min(Xsyn_final)
#mse2=np.std(Xsyn_final)
#snr2=20*np.log10((ima2-imi2)/mse2)
#fprintf(1, 'SNR of SR Image: %f\n',snr3)
#fprintf(1, 'Number of Pixels: %f\n',numel(im2bw(Xsyn3(:))))
#fprintf(1, 'Number of True Pixels: %f\n',sum(im2bw(Xsyn3(:))))
#fprintf(1, 'Maximum Pixel Value relative to Orginial: %f\n',max(Xsyn3(:))) % Max Pixel value
#fprintf(1, 'Reconstruction Mean-square error: %f\n', MSE);
#fprintf(1, 'Reconstruction PSNR: %f\n',PSNR1)