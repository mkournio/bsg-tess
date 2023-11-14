from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.xmatch import XMatch
from astropy.table import Table,  hstack, vstack
import numpy as np
from astropy.io import ascii
from constants import *

#Functions for tables and arrays 

def float_list(a):
	return [float(x) for x in a]

def to_deg(ra,dec):
	
	c = SkyCoord(ra,dec, unit=(u.hourangle, u.deg))

	return c.ra.degree, c.dec.degree

def mask_outliers(m_arr, m = 3.):

    x = np.ma.copy(m_arr)	

    d = np.abs(x - np.ma.median(x))
    mdev = np.ma.median(d)
    s = d / (float(mdev) if mdev else 1.)

    x.mask[s>m] = True

    return x
	

def fill_array(a, fillvalue = 0):

	import itertools

	return np.array(list(itertools.izip_longest(*a, fillvalue=fillvalue))).T

def XmExtCol(ra,dec,ext_files,ext_col,fill_nan_viz=None,**kwargs):

	coord = Table({'RA':ra,'DEC':dec})

	ext_tabs = []
	tflag = len(ext_files)
	for t in reversed(ext_files):
		ext_tab = ascii.read(t)
		ext_tab['f_'+ext_col] = tflag * np.ones(len(ext_tab))
		ext_tab['RA'], ext_tab['DEC'] = to_deg(ext_tab['RA'],ext_tab['DEC'])
		ext_tabs.append(ext_tab)
		tflag -= 1

	ext_glob = vstack(ext_tabs)

	if not ext_col in ext_glob.columns:
			raise Exception('Column name in %s does not exist or is not provided' % ext_file)

	if fill_nan_viz is not None :
		viz_query = XMatch.query(cat1=coord, cat2=fill_nan_viz,  max_distance=2*u.arcsec, \
				     colRA1='RA', colDec1='DEC')
 		viz_col = kwargs.get('viz_col')
		if not viz_col in viz_query.columns:
			raise Exception('Column name in queried Vizier catalogue does not exist or is not provided')
		viz_query = viz_query[['RA','DEC','angDist',viz_col]]
		viz_tab = hstack([coord,match_tabs(coord,viz_query)])

	xmcol = Table(names=[ext_col,'f_'+ext_col],dtype=('f8', 'i2')) 
	for line1 in coord:		
		e_matched = False
		v_matched = False
		match = [0,0]
		for line2 in ext_glob:
			if (line1['RA']==line2['RA']) and (line1['DEC']==line2['DEC']) and \
					not np.ma.is_masked(line2[ext_col]) :
				e_matched = True
				match = [line2[ext_col],line2['f_'+ext_col]]
		if e_matched: 
			xmcol.add_row(match)
		elif fill_nan_viz is not None:
			for line_viz in viz_tab:
				if (line1['RA']==line_viz['RA']) and (line1['DEC']==line_viz['DEC']) and \
					not np.ma.is_masked(line_viz[viz_col]) :
					v_matched = True
					match = [line_viz[viz_col]]
			if v_matched: xmcol.add_row(match)

		if not e_matched and not v_matched : 
			try:
				xmcol.add_row([float(kwargs.get('replace_nan')),0])
			except:
				xmcol.add_row(match,mask=[True,True])

	return xmcol
