;----------------------------------------------------------------------------------------------------------------------------------------------------------
; Spatially Varying Model (SVM)
;----------------------------------------------------------------------------------------------------------------------------------------------------------
;
; This function generates a probability density function using
; Monte-Carlo-Markov-Chains technique, from the target itself.
; Best for using for the artificial star experiment, the PDF 
; provides by this function can be used for the location of 
; the artificial stars in your imaging data.
; INPUTS:
;   idim,jdim: (integer) dimensions of your image
;   xarrcat,yarrcat,magarrcat: (3 arrays with same lengths) position and magnitudes of the detected sources in observed catalogue
;   ArtificialMag: (float value) the magnitude of the artificial star which is going to be added on the original image.
;   mcmcPAR: [N,Sigma] parameters for the MCMC calculation, N is the number of generated positions
;       suggesting value is above 10^5, Sigma is the sigma of the gaussian function using in the MCMC technique to 
;       search for the random positions (suggesting value is 20 for the analysis of the image with angular resolution (FWHM of the PSF) of about 2.5 pixels
;       Note that, if the imaging field is not crowded, then the larger value of sigma is needed, so that MCMC technique could search
;       for positions in a larger domain.
; OUTPUT: is the two dimension array providing the PDF of the positions for artificial stars (X,Y). the length of X,Y positions is N in mcmcPAR      
; 
;  User can randomly sample the positions of the artificial star with magnitude of ArtificialMag among the provided X,Y
;  
;  Calling sequence:
;  
;   myoutput=SVM(idim,jdim, xarrcat, yarrcat, magarrcat, ArtificialMag,[N,Sigma]) 
;   Xoutput=myoutput[*,0]
;   Youtput=myoutput[*,1]
;
; Example:
; observing image size=800x800
; position and magnitude of the detected sources in the catalogue: Xarr,Yarr,mag
; artificial star magnitude = 17.0
; we want to generate 10^6 positions
;
;  myoutput=SVM(800,800,Xarr,Yarr,mag,17.0,[1000000,20])
;  openw,lun,'test',/get_lun
;  for ii=0,(size(myoutput))[1]-1 do printf,lun,myoutput[ii,0],myoutput[ii,1]
;  close,lun
;  free_lun,lun
;
; If you plot 'test' columns1,2, it should be resembeling the shape of your cluster (image) where it is equal or brighter than magnitude 17.
;----------------------------------------------------------------------------------------------------------------------------------------------------------
; Written by Z. Khorrami
; Updated Aug. 2019
; Contact: KhorramiZ@cardiff.ac.uk
;----------------------------------------------------------------------------------------------------------------------------------------------------------


FUNCTION SVM, idim, jdim, xarrcat, yarrcat, magarrcat, ArtificialMag, mcmcPAR

  narrcat=(size(xarrcat))[1]
  if ( (size(yarrcat))[1] ne narrcat and (size(magarrcat))[1] ne narrcat) then begin
    print, 'X,Y,magnitude array should have the same dimention!!!! '
    goto, endoffunction
  endif

  Nnearest=5
  Nmcmc=mcmcPAR[0] ;10000
  Sigmamcmc=mcmcPAR[1] ;20

  image4mcmc=fltarr(idim,jdim)

  newxarr=fltarr(narrcat)
  newyarr=fltarr(narrcat)
  newxarr[*]=100.*idim
  newyarr[*]=100.*jdim

  for ii=0,narrcat-1 do if (magarrcat[ii] le ArtificialMag) then begin
    newxarr[ii]=xarrcat[ii]
    newyarr[ii]=yarrcat[ii]
  endif

  for ii=0,idim-1 do for jj=0,jdim-1 do begin
    one4R=(sort(sqrt((ii-newxarr)^2.+(jj-newyarr)^2.)))[Nnearest-1]
    RR=(ii-newxarr[one4R])^2.+(jj-newyarr[one4R])^2.
    image4mcmc[ii,jj]=1./RR
  endfor

  mymcmc1=mymcmc(image4mcmc,Nmcmc,sigmamcmc)
  writefits,'idl16.fits',image4mcmc
  RETURN,mymcmc1
  endoffunction:

END

FUNCTION mymcmc, img, nn, sigma

  XYarr=fltarr(nn,2)

  idim=(size(img))[1]
  jdim=(size(img))[2]

  newimg=fltarr(idim,jdim)
  newimg[*,*]=img[*,*,0]

  ;seed0=13646
  x0=idim/2.;idim*randomu(seed0)
  y0=jdim/2.;jdim*randomu(seed0)
  if (floor(x0) gt idim-1) then x0=x0-(idim-1)
  if (floor(y0) gt jdim-1) then y0=y0-(jdim-1)

  count=0

  for ii=0,nn do begin
    startagain:
    array=randomn(seed,2)
    x1=sigma*array[0]+x0
    y1=sigma*array[1]+y0
    if (floor(x1) gt idim-1) then x1=x1-(idim-1)
    if (floor(y1) gt jdim-1) then y1=y1-(jdim-1)

    if (floor(x1) lt 0) then x1=(idim-1)+x1
    if (floor(y1) lt 0) then y1=(jdim-1)+y1

    if (newimg[floor(x1),floor(y1)] gt newimg[floor(x0),floor(y0)]) then goto,acceptpoint else begin
      random1=randomu(seed)
      if (newimg[floor(x1),floor(y1)]/newimg[floor(x0),floor(y0)] ge random1) then goto,acceptpoint else goto,startagain
    endelse

    acceptpoint:
    XYarr[count,0] = x1
    XYarr[count,1] = y1
    count=count+1
    x0=x1
    y0=y1
    if (count ge nn) then goto,exitloop
  endfor
  exitloop:

  RETURN,XYarr
END
  
