from astropy import units as u
from astroquery.xmatch import XMatch
import pandas as pd
from astropy.table import Table, hstack
import numpy as np
from constants import *
from scipy.interpolate import interp1d
import os
from astropy.io import ascii

# SEDs, photometry, conversion

def sed_read(temp,logg,models,wave):

	# Read SED and interpolate flux at given wavelengths. 


	if models == 'kurucz':
		sed_name = 't%05dg%02d.flx' % (temp,logg)
  		tab = np.genfromtxt(KURUCZ_SED_PATH+sed_name, names = "wave, flux")
	elif models == 'powr':
		sed_name = 'ob-i_%2d-%2d_sed.txt' % (1e-3*temp,logg)
  		tab = np.genfromtxt(POWR_SED_PATH+sed_name, names = "wave, flux")
		tab['wave'] = 10 ** tab['wave']
		# Conversion into surface flux - raw is flux received at D = 10 pc
		for l in open(POWR_SED_PATH+'modelparameters.txt'):
			if '%2d-%2d'% (1e-3*temp,logg) in l :
				lstar = 10**float(l.split()[5])
				rad2 = lstar * ((float(l.split()[1])/5778.)**-4)
				scalar = rad2 * (1e-1 *  RSUN_TO_PC)**2
		tab['flux'] = (10 **  tab['flux']) * (scalar**-1)
		
	f = interp1d(tab['wave'],tab['flux'],kind='linear')

 	return f(wave)

def getPhot(ra,dec,photo_cats):

	# Retrieve photometry from directory catalogs
	# Data should be a Table containing columns named RA and DEC (in degrees)
	# Return photometric table (in mag) and filter tab dictionary

	print 'Generating photometry from queried catalogues'

	coord = Table({'RA':ra,'DEC':dec})

	filt_tab = {}
	integr=[coord]; table_names=['INPUT']
	for k, v in photo_cats.items():
		query = XMatch.query(cat1=coord, cat2='vizier:'+v['cat'],  max_distance=v['rad']*u.arcsec, \
				     colRA1='RA', colDec1='DEC')
		#print query.columns
		filt = v['flt']
		query = query[['RA','DEC','angDist']+filt]

		if k == 'Merm91': 
			query['Bmag'] = query['Vmag'] + query['B-V']
			query['e_Bmag'] = np.sqrt((query['e_Vmag']**2) + (query['e_B-V']**2))
			query['Umag'] = query['Bmag'] + query['U-B']
			query['e_Umag'] = np.sqrt((query['e_Bmag']**2) + (query['e_U-B']**2))		
			filt = ['Umag','Bmag','Vmag','e_Vmag','e_Bmag','e_Umag']
	
		filt_tab[k] = {}
		filt_tab[k]['f'] = [c for c in filt if not c.startswith('e_')]
		filt_tab[k]['e'] = ['e_'+c for c in filt_tab[k]['f']]
		filt_tab[k]['l'] = [FLZ_DIR[c]['l'] for c in filt_tab[k]['f']]
		filt_tab[k]['z'] = [FLZ_DIR[c]['z'] for c in filt_tab[k]['f']]
		filt_tab[k]['mrk'],filt_tab[k]['fit'] = v['mrk'],v['fit']

		integr.append(match_tabs(coord,query))
		table_names.append(k)

	photo_tab = hstack(integr,uniq_col_name='{col_name}/{table_name}',table_names=table_names)	
	ascii.write(photo_tab, 'photo.dat', overwrite=True)

	return photo_tab, filt_tab

def match_tabs(tab1,tab2):

	merged = Table(tab2[0:0])
	for line1 in tab1:
		min_ang = 99
		match = np.ones(len(tab2.columns))
		for line2 in tab2:
			if (line1['RA']==line2['RA']) and (line1['DEC']==line2['DEC']) and (line2['angDist'] < min_ang):
				min_ang = line2['angDist']
				match = line2
		if min_ang == 99:
			merged.add_row(match,mask=np.ones(len(tab2.columns)))
		else:
			mask = [np.ma.is_masked(match[c]) for c in tab2.columns]
			merged.add_row(match,mask=mask)
				
	del merged['RA','DEC','angDist']

	return merged

def modelTab(filt_dict,models='kurucz'):

	#Creates tab with synthetic fluxes from LTE ATLAS9 or NLTE PoWR grids
	#It is generated over a filter dictionary containing effective wavelengths of bands

	print 'Generating synthetic photometry with %s models' % models

	synth_l = []
	for v in filt_dict.values(): synth_l += v['l']
	synth_l = sorted(set(synth_l))

	mod_tab = Table(names=['t','g']+synth_l)

	if models == 'kurucz':
		grid_teffs = sorted(set([float(f[1:6]) for f in os.listdir(KURUCZ_SED_PATH) if f.endswith('flx')]))
		grid_loggs = sorted(set([float(f[7:9]) for f in os.listdir(KURUCZ_SED_PATH) if f.endswith('flx')]))
	elif models == 'powr':
		grid_teffs = sorted(set([1e+3*float(f[5:7]) for f in os.listdir(POWR_SED_PATH) if f.startswith('ob')]))
		grid_loggs = sorted(set([float(f[8:10]) for f in os.listdir(POWR_SED_PATH) if f.startswith('ob')]))

	for t in grid_teffs:
		for g in grid_loggs:
			try:
					sed = sed_read(t,g,models,synth_l)
					mod_tab.add_row(np.append([t,g],sed))				
			except:
				pass

	return mod_tab

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
