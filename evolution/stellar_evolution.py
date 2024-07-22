from plot_methods import PanelTemplate, colorbar
from astropy.table import Table
from matplotlib.collections import LineCollection
from constants import *
import numpy as np
import os
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier
import pickle
from scipy.interpolate import griddata

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments

class EvolutionaryModels(object):

	def __init__(self, group = 'Geneva', gal = 'MW', rotation = True, load_tr_pickle = False):

		if load_tr_pickle :

			self.tr_tab = pickle.load(open(PICKLE_PATH+'tr.pkl','rb'))
			print 'Loaded pickle: evolutionary models'

		else :
			self.validate()
			catalog = EVOLUTION[group][gal]
			query =	Vizier(columns=['**'], row_limit = -1).get_catalogs(catalog=catalog)
			
			self.tr_tab = query[0]
			if rotation : 
				self.tr_tab = self.tr_tab[self.tr_tab['Rot'] == 'r']

			pickle.dump(self.tr_tab,open(PICKLE_PATH+'tr.pkl','wb'))
			print 'Saved pickle: evolutionary models'

		if group == 'Geneva' :
			self.tr_tab['logg'] = 4 * self.tr_tab['logTe'] - np.log10(self.tr_tab['Gedd'] * Ccm / (KAPPA * SBOLTZ))



		return

	def validate(self):
		pass

	def interp2d(self, key1, key2, key_interp, x, y, post_rsg = False, **kwargs):

		shr_tab = self.tr_tab.copy()
		if not post_rsg:
			mask = shr_tab['Time'] < np.array([PRSG_AGE.get('%s' % int(j), 1e+13) for j in shr_tab['Mini'] ])
			shr_tab = shr_tab[mask]	

		##print shr_tab[key1,key2,key_interp].pprint(max_lines=-1)

		return griddata((shr_tab[key1],shr_tab[key2]),shr_tab[key_interp],(x,y),**kwargs)


	def plot_spectroHR(self, data, xkey = 'TEFF', ykey = 'LOGG', post_rsg = False, hold = False, **kwargs):	

		pl_shr = PanelTemplate(inv_x = True, inv_y = True, x_label = xkey, y_label = ykey, **kwargs)
		ax_shr = pl_shr.PanelAx()

		shr_tab = self.tr_tab.copy()
		if not post_rsg:
			mask = shr_tab['Time'] < np.array([PRSG_AGE.get('%s' % int(x), 1e+12) for x in shr_tab['Mini'] ])
			shr_tab = shr_tab[mask]	
		for l in shr_tab:
			if l['Line'] == 1 and l['Mini'] > 2.4:
				ax_shr.text(l['logTe']+0.01, l['logg']+0.1, int(l['Mini']), size = 10, rotation=90, weight = 'bold')

		flag_fraser = data['f_'+xkey] == 3
		flag_haucke = data['f_'+xkey] == 2
		flag_mattias = data['f_'+xkey] == 1
			
		shr_tab = self._mask_tab(shr_tab, 'Line', 1)
		ax_shr.plot(shr_tab['logTe'],shr_tab['logg'],'k')

		ax_shr.plot(data[xkey][flag_fraser],data[ykey][flag_fraser],'g^', ms = SHR_MS)
		ax_shr.plot(data[xkey][flag_haucke],data[ykey][flag_haucke],'bs', ms = SHR_MS)
		ax_shr.plot(data[xkey][flag_mattias],data[ykey][flag_mattias],'co', ms = SHR_MS)

		#axis.errorbar(tab[k2][flag_haucke],tab[k1][flag_haucke],xerr=tab['e_'+k2][flag_haucke],ecolor='r',**e_kwargs)	
		#axis.errorbar(tab[k2][flag_fraser],tab[k1][flag_fraser],xerr=tab['e_'+k2][flag_fraser],ecolor='g',**e_kwargs)
		
		ax_shr.set_ylim(SHR_YLIM)
		ax_shr.set_xlim(SHR_XLIM)

		if not hold :
			pl_shr.PanelClose()

			return
		else: 
			self.ax = ax_shr
			self.panel = pl_shr

			return 


	def plot(self, key1, key2, key3, **kwargs):

		tr_tab = self.tr_tab.copy()

		mask = tr_tab['Time'] < np.array([PRSG_AGE.get('%s' % int(x), 1e+12) for x in tr_tab['Mini'] ])
		tr_tab = tr_tab[mask]

		fig, ax = plt.subplots()
		cntr = plt.tricontourf(tr_tab[key1],tr_tab[key2],tr_tab[key3], **kwargs)
		fig.colorbar(cntr)
		plt.gca().invert_xaxis()

		plt.show()

		return

	def _mask_tab(self, tab, key, val):

		mtab = tab.copy()
		r = 0
		while r < len(mtab):
			if mtab[r][key] == val :
				for c in mtab.columns :
					mtab[r][c] = np.ma.masked
			r += 1

		return mtab


class getTracks(object):

	def __init__(self, prsg = True):

		self.prsg = prsg

		return

	@property
	def SMC(self):
		return self._read('SMC')

	@property
	def MW(self):
		return self._read('MW')

	def _read(self,gal):

		dict_tr = {}
		_tracks = [x for x in os.listdir(EVOL_TRACKS_PATH) if gal in x]

		for tr in _tracks:
			m_ini = tr.split('_')[0][1:]
			n, m_0, time, mass, logl, logt, vcrit, veq, edd = np.loadtxt(EVOL_TRACKS_PATH + tr, comments='#', unpack = True)

			if not self.prsg:
				mask = time < PRSG_AGE.get(m_ini, 1e+12)
				time = time[mask]; mass = mass[mask]; logl = logl[mask]; logt = logt[mask]; vcrit = vcrit[mask]; veq = veq[mask]

			dict_tr[m_ini] = Table({'TIME' : time, 'M_ACT' : mass, 'LOGL' : logl, 'LOGT' : logt, 'VCRIT' : vcrit, 
						'V_EQ' : veq}, dtype=('f8','f8','f8','f8','f8','f8')) 

  		return dict_tr 	

		
class HRdiagram(PanelTemplate):

	def __init__(self, indata = None, tkey = 'TEFF', cbar = False, lkey = 'LOGL', 
				hold = False, **kwargs):


		self.tkey = tkey
		self.lkey = lkey
		self.cbar = cbar
		super(HRdiagram,self).__init__(inv_x = True, x_label = tkey, y_label = lkey, params = PLOT_PARAMS['panel'], **kwargs)
		self.ax = self.PanelAx()
		self.ax.set_ylim(HR_YLIM)
		self.ax.set_xlim(HR_XLIM)

		if isinstance(indata,Table):

			self.indata = indata
			self._validate_hr()
			self._plot_hr()

		self._plot_tracks()

		if not hold : 
			self.PanelClose()
	
	def _validate_hr(self):

		if not all([x in self.indata.columns for x in [self.tkey,self.lkey]]):
			raise KeyError("Temperature and/or luminosity keys, %s and %s, not found" % (self.tkey,self.lkey))
		return

	def _plot_hr(self):
		
		if not self.cbar :
			colors = 'k'
		else:
			self.cmap, self.norm = colorbar(self.fig, 
				vmin=np.nanmin(indata[self.cbar]), vmax=np.nanmax(indata[self.cbar]), label=STY_LB[cbar])
			colors = np.array([self.cmap(self.norm(k)) for k in self.indata[self.cbar]])


		mask_acyg = (self.indata['VART'] == 'ACYG') | (self.indata['VART'] == 'EB+ACYG')| (self.indata['VART'] == 'ACYG/GCAS')
		mask_sdor = self.indata['VART'] == 'SDOR/L'
		mask_ot = (self.indata['OTTEFF'] == 1) | (self.indata['OTS_LOGL'] == 1)
		
		p_ot, = self.ax.plot(self.indata[self.tkey][mask_ot], self.indata[self.lkey][mask_ot], 'c*', ms = HR_MS+2)
		self.ax.plot(self.indata[self.tkey][~mask_ot], self.indata[self.lkey][~mask_ot], 'bo', ms = HR_MS-3)
		p_acyg, = self.ax.plot(self.indata[self.tkey][mask_acyg], self.indata[self.lkey][mask_acyg], 'bo', mew = 1.5, ms = HR_MS+5, mfc = 'none')
		p_sdor, = self.ax.plot(self.indata[self.tkey][mask_sdor], self.indata[self.lkey][mask_sdor], 'gs', mew = 1.5, ms = HR_MS+5, mfc = 'none')
		self.ax.legend([p_sdor,p_acyg,p_ot],[r's Dor variable', r'$\alpha$ Cygni variables', r'$f_{i}$ > Q3 + 1.5 $\times$ IQR'],labelspacing=.9)	

		return

	def _plot_tracks(self):

		# HD LIMIT
		self.ax.plot(HD_LMT['TEFF'],HD_LMT['LOGL'], lw = HD_LMT['lw'], color = HD_LMT['c'])

		# EVOL TRACKS
		MW_tr = getTracks().MW
		for m in MW_tr:
			self.ax.plot(MW_tr[m]['LOGT'],MW_tr[m]['LOGL'],'k',linewidth=0.6)
			if int(m) > 10: self.ax.text(MW_tr[m]['LOGT'][0]+0.05,MW_tr[m]['LOGL'][0]-0.02,'%s M$_{\odot}$' % int(m), weight= 'bold')

		return
