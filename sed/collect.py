from astropy.table import Table, hstack
from astroquery.xmatch import XMatch
from astropy import units as u
from astropy.io import ascii
from constants import *
from functions import *
from tables import *
import numpy as np
import os

class PhotoCollector(object):

	# Retrieve photometry from directory catalogs
	# Data should be a Table containing columns named RA and DEC (in degrees)
	# Attributes returned; photometric table (in mag) and filter dictionary.

	def __init__(self, input_tab):

		self.filt_dict = {}
		self.input_tab = input_tab
		self.photo_cats = {k:v for k, v in PHOTO_CATS.items() if v['act'] > 0}
		self._collect_phot()

		return

	def _collect_phot(self):

		coord = Table({'RA':self.input_tab['RA'],'DEC':self.input_tab['DEC']})

		integr=[self.input_tab]; table_names=['INPUT']
		for k, v in self.photo_cats.items():

			query = XMatch.query(cat1=coord, cat2='vizier:'+v['cat'],  max_distance=v['rad']*u.arcsec, \
					     colRA1='RA', colDec1='DEC')
			filt = v['flt']
			query = query[['RA','DEC','angDist']+filt]

			if k == 'Merm91': 
				query['Bmag'] = query['Vmag'] + query['B-V']
				query['e_Bmag'] = np.sqrt((query['e_Vmag']**2) + (query['e_B-V']**2))
				query['Umag'] = query['Bmag'] + query['U-B']
				query['e_Umag'] = np.sqrt((query['e_Bmag']**2) + (query['e_U-B']**2))		
				filt = ['Umag','Bmag','Vmag','e_Vmag','e_Bmag','e_Umag']
	
			self.filt_dict[k] = {}
			self.filt_dict[k]['f'] = [c for c in filt if not c.startswith('e_')]
			self.filt_dict[k]['e'] = ['e_'+c for c in self.filt_dict[k]['f']]
			self.filt_dict[k]['l'] = [FLZ_DIR[c]['l'] for c in self.filt_dict[k]['f']]
			self.filt_dict[k]['z'] = [FLZ_DIR[c]['z'] for c in self.filt_dict[k]['f']]
			self.filt_dict[k]['mrk'],self.filt_dict[k]['fit'] = v['mrk'],v['fit']

			integr.append(coord_near_matches(coord,query))
			table_names.append(k)

		self.photo_tab = hstack(integr,uniq_col_name='{col_name}/{table_name}',table_names=table_names)	

		return

	def save_tab(self):

		ascii.write(self.photo_tab, 'phototab.dat', overwrite=True)

		return


