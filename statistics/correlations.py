from plot_methods import GridTemplate
from functions import *
import numpy as np


class Correlation(GridTemplate):

	def __init__(self, ext_x, ext_y, match_keys, **kwargs):

		self.ext_x = ext_x
		self.ext_y = ext_y
		self.match_keys = match_keys

		self.cols_x = [c for c in self.ext_x.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]
		self.cols_y = [c for c in self.ext_y.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]

		super(Correlation,self).__init__(rows_page = len(self.cols_x), cols_page = len(self.cols_y),
						 row_labels = self.cols_x, col_labels = self.cols_y, **kwargs)

		self._corr_panel()
		self.GridClose()	

	def _corr_panel(self):

		from astropy.table import join
		jtab = join(self.ext_x, self.ext_y, keys=self.match_keys)

		try:
			jtab['TEFF'] = np.log10(jtab['TEFF'])
		except:
			pass

		for c1 in self.cols_x :
			for c2 in self.cols_y :

				ax = self.GridAx()

				flag_fraser = jtab['f_'+c2] == 2
				ax.plot(jtab[c2][flag_fraser],jtab[c1][flag_fraser],'g^')

				flag_haucke = jtab['f_'+c2] == 1
				ax.plot(jtab[c2][flag_haucke],jtab[c1][flag_haucke],'ro')

				#ax.plot(jtab[c2],jtab[c1],'k.')				

