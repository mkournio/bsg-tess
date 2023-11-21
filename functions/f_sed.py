

import pandas as pd
from astropy.table import Table, hstack
import numpy as np
from constants import *

import os
from astropy.io import ascii



def mg2flux(mag,emag,zpf,wv,cat):    
 	# Converts magnitudes to fluxes
	# Output in Janskys. 
	# For output in W/m^2/micron append the term  (10**-26) * 2.998E+14 * (wv**-2)
	# For output in erg/s/cm2/A append the term (2.99792458E-05) * (wv**-2)

	mag = pd.to_numeric(np.array(mag), errors = 'coerce')
	emag = pd.to_numeric(np.array(emag), errors = 'coerce')
	zpf = np.array(zpf).astype(np.float)
	wv = np.array(wv).astype(np.float)

  	fl = zpf * 10**(-0.4 * mag) * 2.99792458E-05 * (wv**-2)
	efl = fl * emag * 0.4 * np.log(10)

	if cat == 'AKARI':	#AKARI/MSX give directly Janskys
		fl  = mag *  2.99792458E-05 * (wv**-2)
		efl = emag *  2.99792458E-05 * (wv**-2)
	elif cat == 'MSX6C':
		fl  = mag *  2.99792458E-05 * (wv**-2)
		efl = (2*efl/100.) * fl	
	elif cat == 'IRAS':
		fl  = mag *  2.99792458E-05 * (wv**-2)
		efl = (2*efl/100.) * fl						 

	return fl, efl

def red_coeff(wave,Av,Rv):       
	# Reddening law - cardelli 1998
	# wave should be in microns
	
	Rv = float(Rv)
	x = 1. / wave 
	coeff = np.zeros(len(x))
 
	for i in range(len(x)) :

  		if x[i] < 0.3 :
    			a = 0
    			b = 0

  		if 0.3 <= x[i] < 1.1 :         # Infrared
    			a =  0.574 * (x[i]**1.61)
   			b = -0.527 * (x[i]**1.61)

  		if 1.1 <= x[i] <= 3.3 :        # Optical / Near-Infrared
    			y = x[i] - 1.82   
   			a = 1 + (0.17699*y) - (0.50447*(y**2)) - (0.02427*(y**3)) + (0.72085*(y**4)) + \
				(0.01979*(y**5)) - (0.77530*(y**6)) + (0.32999*(y**7)) 
    			b = (1.41338*y) + (2.28305*(y**2)) + (1.07233*(y**3)) - (5.38434*(y**4)) - (0.62251*(y**5)) + \
				(5.30260*(y**6)) - (2.09002*(y**7))

  		if 3.3 < x[i] <= 8 :           # Ultraviolet 
   
    			if 5.9 <= x[i] <= 8 :
      				fa = -(0.04473*((x[i]-5.9)**2)) - (0.009779*((x[i]-5.9)**3))
      				fb =  (0.21300*((x[i]-5.9)**2)) + (0.120700*((x[i]-5.9)**3))
    			if x[i] < 5.9:
      				fa=0
      				fb=0

    			a =  1.752 - (0.316*x[i]) - (0.104 / float( ((x[i]-4.67)**2) + 0.341)) + fa
    			b = -3.090 + (1.825*x[i]) + (1.206 / float( ((x[i]-4.62)**2) + 0.263)) + fb

  		if 8 < x[i] <= 10 :            # Far UV

    			a = -1.073 - (0.628*(x[i]-8)) + (0.137*((x[i]-8)**2)) - (0.070*((x[i]-8)**3))	
    			b = 13.670 + (4.257*(x[i]-8)) - (0.420*((x[i]-8)**2)) + (0.374*((x[i]-8)**3))	

		coeff[i] = Av * (a + b / Rv)
 
  	return coeff

def planck(w,temp,beta=1.5,**kwargs):

	w = np.array(w,dtype=float)

	#Convert A to m
	w = w * 1e-10

	C1 = 2 * PI * hP * (C**2)
	C2 = hP * C / kB

	B = C1 * (w**-5) * (np.exp(C2 * (w * temp)**-1) - 1)**-1

	if kwargs['kind'] == 'grey' : B *= (w**-beta)

	#Convert W/(m^2 m) to erg/(cm^2 A s)
	return B * 1e-7 
