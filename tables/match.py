from astropy import units as u
from astroquery.xmatch import XMatch
from astropy.table import Table,  hstack, vstack
import numpy as np
from astropy.io import ascii
from constants import *
from functions import *
import pickle

def coord_near_matches(tab1,tab2):

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


class XMatching(object):

	def __init__(self, data_tab, xtabs, xcols, vizier = {}, xrad = 2, load_pickle = False):

		if load_pickle:
			self.matched = pickle.load(open(PICKLE_PATH+'input.pkl','rb'))
			print 'Loaded pickle: input'
		else:
			self.matched = data_tab
			self.coord = Table({'RA': data_tab['RA'] ,'DEC' : data_tab['DEC']})
			self.xrad = xrad
			self.match(xtabs = xtabs, xcols = xcols)

			for k in vizier:
				self.match(xtabs = [str(k)], xcols = [vizier[k]['xcols']], viz_col = vizier[k]['viz_col'])

			if 'TIC' in self.matched.columns: self.matched['TIC'].format = '%10d'
			if 'HDIST' in self.matched.columns:
			 self.matched['HDIST'] = 1./(1e-3*self.matched['HDIST'])	
			 self.matched['HDIST'].mask = self.matched['HDIST'] < 0.
	
			pickle.dump(self.matched,open(PICKLE_PATH+'input.pkl','wb'))


	def match(self, xtabs, xcols, **kwargs):

		for xcol in xcols:

			self.xcol = xcol
			viztype = ['vizier:' in x for x in xtabs]

			if len(xtabs) == 0 :
				raise Exception("External table(s) list is not provided - quitting")
			elif not any(viztype):
				xmcol = self._match_ext_tab(xtabs)
			elif all(viztype):
				xmcol = self._match_vizier_tab(xtabs,**kwargs)
			else :
				raise Exception("Unsupported yet type list")

			self.matched = hstack([self.matched,xmcol])

		return


	def _match_ext_tab(self, xtab_files):

		xtabs = []	
		tflag = len(xtab_files)
		for xtab in reversed(xtab_files):
			xtab = ascii.read(xtab)
			xtab['f_'+ self.xcol] = tflag * np.ones(len(xtab))
			xtab['RA'], xtab['DEC'] = to_deg(xtab['RA'],xtab['DEC'])
			xtabs.append(xtab)
			tflag -= 1

		stacked_xtab = vstack(xtabs)
		if self.xcol not in stacked_xtab.columns:
			raise Exception("Requested column %s not present in the external tables" % self.xcol)

		return self._matching(stacked_xtab)


	def _match_vizier_tab(self,xtabs,**kwargs):

		xcat = xtabs[0]
		viz_query = XMatch.query(cat1 = self.coord, cat2 = xcat,  max_distance = self.xrad * u.arcsec, \
				         colRA1 = 'RA', colDec1 = 'DEC')

		if ('viz_col' not in kwargs) or (kwargs['viz_col'] not in viz_query.columns):
			raise Exception("Vizier: you provided wrong or no column name for match")

		viz_col = kwargs['viz_col']
		viz_query = viz_query[['RA','DEC','angDist',viz_col]]
		viz_query.rename_column(viz_col, self.xcol)
		viz_tab = hstack([self.coord,coord_near_matches(self.coord,viz_query)])
		viz_tab['f_'+self.xcol] = -1 * np.ones(len(viz_tab))

		return self._matching(viz_tab,**kwargs)


	def _matching(self,mtab,**kwargs):

		column_format = kwargs.get('column_format','f8')
		xmcol = Table(names=[self.xcol,'f_'+self.xcol],dtype=(column_format, 'i2')) 

		for line1 in self.coord:		
			matched = False
			match = [0,0]
			for line2 in mtab:
				if (line1['RA']==line2['RA']) and (line1['DEC']==line2['DEC']) and \
					not np.ma.is_masked(line2[self.xcol]) :
					matched = True
					match = [line2[self.xcol],line2['f_'+self.xcol]]
			if matched: 
				xmcol.add_row(match)
			else:
				try:
					xmcol.add_row([float(kwargs.get('replace_nan')),0])
				except:
					xmcol.add_row(match,mask=[True,True])
		return xmcol
		


