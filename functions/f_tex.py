from astropy.table import Table
import numpy as np
from astropy.io import ascii
from constants import *

class TexTab(object):

	def __init__(self):		
		pass

	def TabSample(self, data):

		data_tex = Table(data['STAR','RA','DEC', 'SBTYPE', 'TIC', 'RDMAG'], 
		names = ('Star','RA', 'DEC', 'Sp. type', 'TIC', 'G$_{RP,*}$ - G$_{RP,2}$'))
		data_tex['RA'].format = '2.6f'; data_tex['DEC'].format = '2.6f'
		sect_col = []
		for i,j in zip(data['SSPOC'],data['SFFI']):
			sect = []
			if not np.ma.is_masked(i): sect += [int(x) for x in str(i).split(';')]
			if not np.ma.is_masked(j): sect += [int(x) for x in str(j).split(';')]
			sect_col.append(','.join([str(x) for x in sorted(sect)]))
		data_tex['TESS sectors'] = np.array(sect_col)

		TEX_SAMPLE_TAB['units'] = {'RA' : '(deg)', 'DEC' : '(deg)', 'G$_{RP,*}$ - G$_{RP,2}$' : '(mag)'}
		TEX_SAMPLE_TAB['caption'] = r'\label{tab_sample} Sample of the studied BSGs.'
		TEX_SAMPLE_TAB['col_align'] = 'lrrl|ccr'
 
		ascii.write(data_tex, 'tab_sample.tex', format="latex", latexdict=TEX_SAMPLE_TAB)

		return
