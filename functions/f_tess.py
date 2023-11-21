# Functions related to the analysis of the TESS data

from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import numpy as np
from constants import *
from lightkurve.correctors import DesignMatrix,RegressionCorrector
from scipy import signal
from tables import *

class Gaia(object):

    def __init__(self,tpf,cat='Gaia3'):

        self.tpf = tpf
        if cat=='Gaia3' :
            self.cat = 'I/355/gaiadr3'
        elif cat=='Gaia2' :
            self.cat = 'I/345/gaia2'
        else :
            raise ValueError("Choose between Gaia3 and Gaia2.")

	return

    def _query(self):

        DS3 = Vizier(catalog=self.cat,columns=['*','+_r']); DS3.ROW_LIMIT = -1
        DS3_query = DS3.query_region(SkyCoord(ra=self.tpf.ra,dec=self.tpf.dec,unit=(u.deg,u.deg), frame='icrs'),
                                     radius = max(self.tpf.shape[1:]) * TESS_pix_size * u.arcsec)
        DS3 = DS3_query[0].to_pandas()

        return DS3

    def get_prop(self):

        DS3 = self._query()

        RA_pix,DE_pix = self.tpf.wcs.all_world2pix(DS3.RA_ICRS,DS3.DE_ICRS,0.01)
        RA_pix += self.tpf.column ; DE_pix += self.tpf.row
        RA_pix += 0.5 ; DE_pix += 0.5
        
        cond_box = (RA_pix>self.tpf.column) & (RA_pix<self.tpf.column+self.tpf.shape[1]) & \
                          (DE_pix>self.tpf.row) & (DE_pix<self.tpf.row+self.tpf.shape[2])
        DS3 = DS3[cond_box]
        RA_pix = RA_pix[cond_box]
        DE_pix = DE_pix[cond_box]

        BPRP  = DS3['BP-RP'].values
        RP    = DS3['RPmag'].values
	GaID  = DS3['Source'].values

        return  RA_pix, DE_pix, BPRP, RP, GaID

def find_harmon(v,e):

	'''
	Calculates the harmonics/combinations from a vector of frequencies
	'''
	v = np.array(v)
	ident_v = np.zeros(len(v),dtype=np.dtype('U100'))

	f = np.where(v > 2 * e)[0]
	t = np.where(v <= 2 * e)[0]

	np.put(ident_v,t,'(BTH)')

	ratio = v[f[0]]/v[f[1]]
	error = ratio * e * np.sqrt( (1/v[f[0]]**2) + (1/v[f[1]]**2) )	

	if (abs(ratio - round(ratio)) < error) and (round(ratio) == 1) :
			id0 = '(F)'
			id1 = '(U) F%d' % f[0]
	elif abs(ratio - round(ratio)) < error :
			id0 = '(H) %dF%d' % (round(ratio),f[1])
			id1 = '(F)'
	elif abs((1/ratio) - round(1/ratio)) < error * (1/ratio**2) :
			id0 = '(F)'
			id1 = '(H) %dF%d' % (round(1/ratio),f[0])
	else:
			id0 = '(F)'
			id1 = '(F)'

	ident_v[f[0]] = id0
	ident_v[f[1]] = id1

	for k in f[2:] :

		ind_k = np.where(f==k)[0][0]

		minim = 100.
 		oc1 = '0F'; oc2 = '0F'
 		for i in f[:ind_k-1] :
			for j in f[i+1:ind_k] :
					for c1 in np.arange(-5,10,1):
    						for c2 in np.arange(-5,10,1):

				        		if (abs(c1*v[i] + c2*v[j] - v[k]) < np.sqrt(2) * e) and \
								(abs(c1) + abs(c2) < minim): 

								minim = abs(c1) + abs(c2)
								oc1 = '%sF%s' % (c1,i); oc2 = '%sF%s' % (c2,j)

	 	sum_coeff = abs(float(oc1.split('F')[0]))+abs(float(oc2.split('F')[0]))
		oc1 = trans_freq_text(oc1)
		oc2 = trans_freq_text(oc2)

	 	text = '%s %s' % (oc1,oc2)

	 	if sum_coeff == 0:
			text = '(F)'
	 	elif sum_coeff == 1:
			text = '(U) ' + text
	 	elif sum_coeff > 1 and ('' in [oc1,oc2]):
			text = '(H) ' + text
	 	else:
			text = '(C) ' + text

	 	ident_v[k] = text 

	return ident_v

def trans_freq_text(t):

	if float(t.split('F')[0]) == 0. : 
		return ''
	elif float(t.split('F')[0]) == 1. :
		return t[1:]
	else:
		return t

def read_peak(file_name):

	loc_peaks = []
	with open(file_name) as f:
		glob_peaks = [float(x) for x in f.readline().split()[2:]]
		for line in f:
			loc_peaks.append([float(x) for x in line.split()[2:]])

	return np.array(glob_peaks), np.array(loc_peaks)

def group_consec_lcs(v):

	import more_itertools as mit
	v = sorted(v, key = lambda x: sect_from_lc(x))

	return [list(group) for group in mit.consecutive_groups(v, ordering = lambda x: sect_from_lc(x))]

def sect_from_lc(v):

	return int(v.split('_')[1])

def lc_read(lc_file_path, remove_nans = True):

	time, flux = np.loadtxt(lc_file_path, delimiter=' ', usecols=(0, 1), unpack=True)
	if remove_nans: flux = set_nans(time, flux)

	return time, flux

def read_prg(file_name):

	x, y = np.loadtxt(file_name, delimiter=' ', usecols=(0, 1), unpack=True)

	return x, y
	

def getmask(tpf, star_row, star_col, thres = 0.1):

	mask = np.zeros(tpf[0].shape[1:], dtype='bool')
	mask[star_row][star_col] = True

	flux_matr = tpf[0].flux[0]
	cen_flux = flux_matr[star_row][star_col]

	rad = 1
	mask_size = len(mask)

	for c in np.arange(star_col,mask_size,1) :

			col_break = -1 

			for r in np.arange(star_row,-1,-1):

				max_flx = -1

				if flux_matr[r][c] > thres * cen_flux : 
					mask[r][c] = True
					max_flx = 1
					col_break = 1

				if max_flx < 0. : break

			for r in np.arange(star_row,mask_size,1):

				max_flx = -1

				if flux_matr[r][c] > thres * cen_flux : 
					mask[r][c] = True
					max_flx = 1
					col_break = 1

				if max_flx < 0. : break	

			if col_break < 0. : break	

	for c in np.arange(star_col,-1,-1) :

			col_break = -1 

			for r in np.arange(star_row,-1,-1):

				max_flx = -1

				if flux_matr[r][c] > thres * cen_flux : 
					mask[r][c] = True
					max_flx = 1
					col_break = 1

				if max_flx < 0. : break

			for r in np.arange(star_row,mask_size,1):

				max_flx = -1

				if flux_matr[r][c] > thres * cen_flux : 
					mask[r][c] = True
					max_flx = 1
					col_break = 1

				if max_flx < 0. : break	

			if col_break < 0. : break				
	
	return mask	

def save_mask(mask,filename):

	np.savetxt(filename, mask, fmt="%5i")

	return

def read_mask(filename):

	return np.loadtxt(filename,dtype=bool)

def dmatr(matr, pca_num):

	return DesignMatrix(matr,name='regressors').pca(pca_num).append_constant()

def lccor(tpf, mask, bk_mask, pca_num): 

	rgr = tpf.flux[:, bk_mask]

	lc    = tpf.to_lightcurve(aperture_mask=mask).remove_nans()
	lc    = lc[lc.flux_err > 0]

	dm = dmatr(rgr,pca_num)

	return RegressionCorrector(lc).correct(dm)

def set_nans(time,flux):

	time_dif = time[1:] - time[:-1]; ind_nan = np.argmax(time_dif)
	flux[ind_nan-1:ind_nan+1] = np.nan;
	flux[0] = np.nan; flux[-1] = np.nan;

	return flux

def fit_pol(x,y,yerr,deg,unit = 'mag'):

	idx = np.isfinite(x) & np.isfinite(y)
	x = x[idx]; y = y[idx]; yerr = yerr[idx]

	p = np.poly1d(np.polyfit(x, y, deg))
	yfit = p(x)

	yfit_norm = y/yfit
	e_yfit_norm = yerr/yfit

	dm = -2.5 * np.log10(yfit_norm)
	e_dm = 2.5 * yerr / (LN10*y)

	if unit == 'mag':
		return x, yfit, dm, e_dm
	elif unit == 'norm':
		return x, yfit, ynorm, e_yfit_norm

def fmin_fit(y,deg):
		
	x = np.arange(0,len(y),1)
	f = np.poly1d(np.polyfit(x,y,deg))

	return np.argmin(f(x))

def save_two_col(x,y,filename):


	with open(filename,'w') as f:
		for i,j in zip(x,y):
    			f.write("%s %s\n" % (i,j))
	f.close()

	return

def save_three_col(x,y,z,filename):

	with open(filename,'w') as f:
		for i,j,k in zip(x,y,z):
    			f.write("%s %s %s\n" % (i,j,k))
	f.close()

	return

def frequencies(tpf, mask, bk_mask, pca = 1, peak_thres = 3):

	print '..global peaks..'
	low_freq = 2 / (tpf.time[-1] - tpf.time[0])
	cor_lc = lccor(tpf,mask,bk_mask,pca).normalize()
	glob_pg = cor_lc.to_periodogram()  
	glob_peaks,_ = signal.find_peaks(glob_pg.power)
	glob_ps = [p for p in glob_peaks if glob_pg.power[p] > (np.nanmedian(glob_pg.power) + peak_thres*np.nanstd(glob_pg.power))]
	glob_ps_freq = glob_pg.frequency[glob_ps].value
	glob_ps_freq = glob_ps_freq[glob_ps_freq>low_freq]
			
	mask_h, mask_l = tpf[0].shape[1:][0], tpf[0].shape[1:][1]
	peaks_array = []
	mask_r = []; mask_c = []
	loc_mask =  np.zeros(tpf[0].shape[1:], dtype='bool')

	print '..local peaks..'
	for r in np.arange(mask_h):
		for c in np.arange(mask_l):

			if mask[r][c] : 

				mask_r.append(r); mask_c.append(c)
				
				loc_mask[r][c] = True

				cor_loc_lc = lccor(tpf,loc_mask,bk_mask,pca).normalize() 
				loc_pg = cor_loc_lc.to_periodogram() 
				loc_peaks,_ = signal.find_peaks(loc_pg.power)
				loc_ps = [p for p in loc_peaks if loc_pg.power[p] > (np.nanmedian(loc_pg.power) + peak_thres*np.nanstd(loc_pg.power))]
				loc_ps_freq = loc_pg.frequency[loc_ps].value
				loc_ps_freq = loc_ps_freq[loc_ps_freq>low_freq]
				peaks_array.append(loc_ps_freq)

				loc_mask[r][c] = False

				#w = np.linspace(0.01, 3, 1000)
    				#pg = signal.lombscargle(cor_loc_lc.time, cor_loc_lc.flux, w, normalize=True)
				#loc_peaks,_ = signal.find_peaks(pg)
				#loc_ps = [p for p in loc_peaks if pg[p] > peak_thres * np.nanmean(pg)]
				#peaks_array.append(w[loc_ps])

	filename = '%s%s_%s_PEAKS%s_pca%s' % (tpf.ra, tpf.dec, tpf.sector, peak_thres,pca)
	with open(filename,'w') as f:
		s = " ".join(map(str, glob_ps_freq))
		f.write("G G %s\n" % s)	
		for r,c,p in zip(mask_r,mask_c,peaks_array) :
			s = " ".join(map(str, p))
			f.write("%s %s %s\n" % (r,c,s))
	f.close()

	return glob_pg, glob_ps_freq, peaks_array


def freqs_to_tex_tab(d_fund, d_harm, d_comb, file_):

	with open(file_, 'w') as tex: 

		tex.write(TEX_FREQ_TAB['PRE'])

		for star in d_fund.keys():

			for index_f, sect in enumerate(d_fund[star]['SECT']) :

				for f, ident, sn in zip(d_fund[star]['FREQ'][index_f],d_fund[star]['IDS'][index_f],d_fund[star]['SNR'][index_f]):
					 tex.write("%s & %s & %s & %.3f & %.1f\\\\ \n" % (star,sect,ident,f,sn))

				index_h = d_harm[star]['SECT'].index(sect)
				for f, ident, sn in zip(d_harm[star]['FREQ'][index_h],d_harm[star]['IDS'][index_f],d_harm[star]['SNR'][index_h]):
					 tex.write("%s & %s & %s & %.3f & %.1f\\\\ \n" % (star,sect,ident,f,sn))

				index_c = d_comb[star]['SECT'].index(sect)
				for f, ident, sn in zip(d_comb[star]['FREQ'][index_c],d_comb[star]['IDS'][index_f],d_comb[star]['SNR'][index_c]):
					 tex.write("%s & %s & %s & %.3f & %.1f\\\\ \n" % (star,sect,ident,f,sn))

				tex.write("\\hline\n")

		tex.write(TEX_FREQ_TAB['POST'])
	tex.close()

	return
	

######################################

def roundup(x):
	import math	

	return int(math.ceil(x / 10.0)) * 10

def sect_unmasked(x):
	
	if not np.ma.is_masked(x):
		return [int(i) for i in str(x).split(';')]
	else:
		return []



