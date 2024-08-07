from astropy.table import Table, Column
import numpy as np
from astropy.io import ascii
from constants import *

class TexTab(object):

	def __init__(self):		
		pass

	def TabSample(self, data):

		data_tex = Table(data['STAR','RA','DEC', 'SBTYPE', 'DIST', 'TIC', 'TMAG', 'RDMAG'], 
		names = ('Star','RA', 'DEC', 'Spectral type','$D$', 'TIC', 'TESS mag','G$_{RP,*}$ - G$_{RP,2}$'))

		data_tex['RA'].format = '2.6f'; data_tex['DEC'].format = '2.6f'; data_tex['TESS mag'].format = '.1f'

		if 'HDIST' in data.columns:
		 for d in range(len(data)):
			if data_tex['$D$'][d] is np.ma.masked : data_tex['$D$'][d] = data['HDIST'][d]
		data_tex['$D$'].format = '%.f'

		sect_col = []
		for i,j in zip(data['SSPOC'],data['SFFI']):
			sect = []
			if not np.ma.is_masked(i): sect += [int(x) for x in str(i).split(';')]
			if not np.ma.is_masked(j): sect += [int(x) for x in str(j).split(';')]
			sect_col.append(','.join([str(x) for x in sorted(sect)]))
		data_tex['TESS sectors'] = np.array(sect_col)

		TEX_SAMPLE_TAB['units'] = {'RA' : '(deg)', 'DEC' : '(deg)', 'TESS mag': '(mag)', '$D$' : '(pc)', 'G$_{RP,*}$ - G$_{RP,2}$' : '(mag)'}
		TEX_SAMPLE_TAB['caption'] = r'\label{tab_sample} Sample of the studied BSGs. Distances for indices 0.1,2,6,10 are taken from Hipparcos.'
		TEX_SAMPLE_TAB['col_align'] = 'lrrlc|cccr'
		TEX_SAMPLE_TAB['preamble'] = r'\small\centering'
 
		ascii.write(data_tex, 'tab_sample.tex', format="latex", latexdict=TEX_SAMPLE_TAB)

		return

	def TabPhoto(self,data):

		opt_columns = ['Umag','Bmag/Merm91','Vmag/Merm91','Bmag/NOMAD05','Vmag/NOMAD05','Rmag','BPmag','Gmag','RPmag']
		ir_columns =['Jmag','Hmag','Kmag','S09','S18','B1','B2','A','C','D','E']

		opt_tex = Table([data['STAR']],names=['Star'])
		ir_tex = Table([data['STAR']],names=['Star'])

		TEX_SAMPLE_TAB['preamble'] = r'\scriptsize'
		TEX_SAMPLE_TAB['caption'] = r'\label{tab_photo} Multi-band photometry of the studied stars. Uncertainties are provided in parentheses when available. Units are in given in mag.'
		TEX_SAMPLE_TAB['tabletype'] = 'table'

		self._new_tab(data,opt_columns,opt_tex)
		self._new_tab(data,ir_columns,ir_tex)

		ascii.write(opt_tex, 'tab_photo_opt.tex', format="latex", latexdict=TEX_SAMPLE_TAB)
		ascii.write(ir_tex, 'tab_photo_ir.tex', format="latex", latexdict=TEX_SAMPLE_TAB)

		return

	def TabCalcProp(self,intab):

		columns = ['STAR','A_V','LOGC','LOGW','LOGR0','TAU','GAMMA', 'SVAR','PSI','SKEW','F_ROT','OTTEFF','OTS_LOGL']
		names = ['Star','$A_{V}$', 'log$C$', 'log$W$','log$R_{0}$',r'$\tau$',r'$\gamma$', r'$\sigma$',r'$\psi^2$',r'skw','R','OS','OL']
		sed_tex = Table(intab[columns],names=names)

		for n in names[1:]:
			sed_tex[n].format = '.2f'	
		sed_tex['$\sigma$'].format = '.3f'
		for n in ['R','OS','OL']:
			sed_tex[n].format = '%d'

		TEX_SAMPLE_TAB['preamble'] = r'\small\centering' 
		TEX_SAMPLE_TAB['caption'] = r'\label{tab:cprop} Calculated parameters.'
		TEX_SAMPLE_TAB['tabletype'] = 'table*'
		TEX_SAMPLE_TAB['col_align'] = 'lcc|cccc|ccc|ccc'

		ascii.write(sed_tex, 'tab_cprop.tex', format="latex", latexdict=TEX_SAMPLE_TAB)

		return	

	def TabSED(self,data):

		columns = ['STAR','A_V','LOGC','S_LOGL']
		names = ['Star','$A_{V}$', 'log($R [R_{\odot}]/D [pc]$)', 'log($L/L_{\odot})_{Gaia}$']
		sed_tex = Table(data[columns],names=names)

		for n in names[1:]:
			sed_tex[n].format = '.2f'

		TEX_SAMPLE_TAB['preamble'] = r'\centering'
		TEX_SAMPLE_TAB['caption'] = r'\label{tab_sed} Calculated parameters from the modeling of the SEDs.'
		TEX_SAMPLE_TAB['tabletype'] = 'table'
		TEX_SAMPLE_TAB['units'] = {'$A_{V}$' : '(mag)'}
		TEX_SAMPLE_TAB['col_align'] = 'lcc|c'

		ascii.write(sed_tex, 'tab_sed.tex', format="latex", latexdict=TEX_SAMPLE_TAB)

		return		


	def _new_tab(self,data,cols,ntab):

		for c in cols:
			if 'e_'+c in data.columns:
				new_col = []
				for i in range(len(data)):
					if not np.ma.is_masked(data['e_'+c][i]):
					 new_col.append('%.2f(%.2f)' % (data[c][i],data['e_'+c][i]))	
					elif not np.ma.is_masked(data[c][i]):
					 new_col.append('%.2f' % data[c][i])
					else:
					 new_col.append('$-$')		
				ntab.add_column(Column(new_col, name=c, dtype='U11'))
			else:
				data[c].fill_value = '$-$'
				ntab.add_column(data[c].filled())

		for c in ntab.columns:
			nc = c.replace('mag','')
			nc = c.replace('BP','G$_{BP}$')
			nc = c.replace('RP','G$_{RP}$')
			nc = nc.replace('/Merm91','(M91)')
			nc = nc.replace('/NOMAD05','(N05)')
			ntab.rename_column(c,nc)

		return 

