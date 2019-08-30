#!/usr/bin/env python
import sys
from math import *
import random
import numpy as np
from scipy.linalg import expm
from astropy.io import fits
import os, sys, shutil
import timeit
import math


#---------------------------------------------------------------------------------------------------
# Spatially Varying Model (SVM)
#---------------------------------------------------------------------------------------------------
#
# This function generates a probability density function using
# Monte-Carlo-Markov-Chains technique, from the target itself.
# Best for using for the artificial star experiment, the PDF 
# provides by this function can be used for the location of 
# the artificial stars in your imaging data.
# INPUTS:
#   idim,jdim: (integer) dimensions of your image
#   xarrcat,yarrcat,magarrcat: (3 arrays with same lengths) position and magnitudes of the detected sources in observed catalogue
#   ArtificialMag: (float value) the magnitude of the artificial star which is going to be added on the original image.
#   mcmcPAR: [N,Sigma] parameters for the MCMC calculation, N is the number of generated positions
#       suggesting value is above 10^5, Sigma is the sigma of the gaussian function using in the MCMC technique to 
#       search for the random positions (suggesting value is 20 for the analysis of the image with angular resolution (FWHM of the PSF) of about 2.5 pixels
#       Note that, if the imaging field is not crowded, then the larger value of sigma is needed, so that MCMC technique could search
#       for positions in a larger domain.
# OUTPUT: is the two dimension array providing the PDF of the positions for artificial stars (X,Y). the length of X,Y positions is N in mcmcPAR      
# 
#  User can randomly sample the positions of the artificial star with magnitude of ArtificialMag among the provided X,Y
#  
#  Calling sequence:
#  
#   myoutput=SVM(idim,jdim, xarrcat, yarrcat, magarrcat, ArtificialMag,[N,Sigma]) 
#   Xoutput=myoutput[*,0]
#   Youtput=myoutput[*,1]
#
# Example:
# observing image size=800x800
# position and magnitude of the detected sources in the catalogue: Xarr,Yarr,mag
# artificial star magnitude = 17.0
# we want to generate 10^6 positions
# lun=open('test',"w")
#  myoutput=SVM(800,800,Xarr,Yarr,mag,17.0,[1000000,2.])
#  for ii in range(Nmcmc): lun.write("%13.4f %13.4f \n" %(myoutput[ii,0],myoutput[ii,1]))
#
# If you plot 'test' columns1,2, it should be resembeling the shape of your cluster (image) where it is equal or brighter than magnitude 17.
#
#---------------------------------------------------------------------------------------------------
# Written by Z. Khorrami
# Updated Aug. 2019
# Contact: KhorramiZ@cardiff.ac.uk
#---------------------------------------------------------------------------------------------------


def mymcmc(img, nn, sigma):

  XYarr=np.full((int(nn),int(2)),0.0)
  jdim,idim=img.shape

  newimg=np.full((int(jdim),int(idim)),0.0) #fltarr(idim,jdim)
  newimg[:,:]=img[:,:]
  #seed0=13646
  x0=idim/2. #idim*randomu(seed0)
  y0=jdim/2. #jdim*randomu(seed0)
  if (math.floor(x0) > idim-1): x0=x0-(idim-1)
  if (math.floor(y0) > jdim-1): y0=y0-(jdim-1)

  count=0

  while (nn > count):   
    array=2.*np.random.rand(2)-1.
    x1=sigma*array[0]+x0
    y1=sigma*array[1]+y0
    if (math.floor(x1) > idim-1): x1=x1-(idim-1)
    if (math.floor(y1) > jdim-1): y1=y1-(jdim-1)

    if (math.floor(x1) < 0): x1=(idim-1)+x1
    if (math.floor(y1) < 0): y1=(jdim-1)+y1

    if (newimg[int(math.floor(y1)),int(math.floor(x1))] > newimg[int(math.floor(y0)),int(math.floor(x0))]): 
      XYarr[count,0] = x1
      XYarr[count,1] = y1
      count=count+1
      x0=x1
      y0=y1
    else: 
      random1=np.random.rand(1)
      if (newimg[int(math.floor(y1)),int(math.floor(x1))]/newimg[int(math.floor(y0)),int(math.floor(x0))] > random1):
        XYarr[count,0] = x1
        XYarr[count,1] = y1
        count=count+1
        x0=x1
        y0=y1
  return XYarr


def SVM(idim, jdim, xarrcat, yarrcat, magarrcat, ArtificialMag, mcmcPAR):

  narrcat=len(xarrcat)
  if ( len(yarrcat) != narrcat and len(magarrcat) != narrcat):
    print 'X,Y,magnitude array should have the same dimention!!!! '
    close()
  

  Nnearest=5
  Nmcmc=mcmcPAR[0] #10000
  Sigmamcmc=mcmcPAR[1] #20

  image4mcmc=np.full((int(jdim),int(idim)),0.0)

  newxarr=np.full(int(narrcat),100.*idim)
  newyarr=np.full(int(narrcat),100.*jdim)

  for ii in range(narrcat):
      if (magarrcat[ii] <= ArtificialMag):
          newxarr[ii]=xarrcat[ii]
          newyarr[ii]=yarrcat[ii]

  for ii in range(idim):
    for jj in range(jdim):
      rdis=np.sqrt((newxarr-np.array(ii))**2.+(newyarr-np.array(jj))**2.)
      one4R=(np.argsort(rdis))[int(Nnearest)-1]
      RR=(ii-newxarr[one4R])**2.+(jj-newyarr[one4R])**2.
      image4mcmc[jj,ii]=1./RR
  hdu = fits.PrimaryHDU(image4mcmc)
  hdu.writeto('testpy.fits')
  mymcmc1=mymcmc(image4mcmc,Nmcmc,Sigmamcmc)

  return mymcmc1

