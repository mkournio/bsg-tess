from plot_methods import PanelTemplate, colorbar
from astropy.table import Table
from matplotlib.collections import LineCollection
from constants import *
import numpy as np
import os

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments


class getTracks(object):

	def __init__(self):
		pass

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

			time, m_act, logl, lteff, veq = [], [], [], [], []
			with open(EVOL_TRACKS_PATH + tr) as file_:
				for line in file_:
				    if not line.startswith('#'):
				      line = line.split()	
				      time.append(line[1])
				      m_act.append(line[2])
				      logl.append(line[3])
				      lteff.append(line[4])
				      veq.append(line[5])				
			file_.close()

			dict_tr[m_ini] = Table({'TIME' : time, 'M_ACT' : m_act, 'LOGL' : logl, 
					       'LTEFF' : lteff, 'V_EQ' : veq}, dtype=('f8','f8','f8','f8','f8')) 

  		return dict_tr 	
		
class HRdiagram(PanelTemplate):

	def __init__(self, data = None, tkey = 'TEFF', cbar_key = 'MASS', lkey = 'LOGL', tlines = [], llines = [], 
				hold = False, **kwargs):


		self.tkey = tkey
		self.lkey = lkey
		self.cbar_key = cbar_key
		self.tlines = tlines
		self.llines = llines
		super(HRdiagram,self).__init__(inv_x = True, x_label = 'LTEFF', y_label = lkey, params = PLOT_PARAMS['panel'], **kwargs)
		self.ax = self.PanelAx()

		self.ax.set_ylim(HR_YLIM)
		self.ax.set_xlim(HR_XLIM)		

		self.cmap, self.norm = colorbar(self.fig, vmin=1, vmax=60, label=r'M [M$_{\odot}$]')
		#self.cmap, self.norm = colorbar(self.fig, vmin=0, vmax=15, label=r't [Myr]')

		self._plot_tracks()

		if isinstance(data,Table):

			self.data = data
			self._validate_hr()
			self.data[self.tkey] = np.log10(self.data[self.tkey])
			self._plot_hr()

		if not hold : 
			self.PanelClose()
	
	def _validate_hr(self):

		if not all([x in self.data.columns for x in [self.tkey,self.lkey]]):
			raise KeyError("Temperature and/or luminosity keys, %s and %s, not found" % (self.tkey,self.lkey))
		return

	def _plot_hr(self):

		msizes = 60 #np.array([x**2 for x in self.data['NFFREQ']])
		colors = np.array([self.cmap(self.norm(k)) for k in self.data[self.cbar_key]])
		self.ax.scatter(self.data[self.tkey], self.data[self.lkey], s = msizes, c = colors, edgecolor='k')		

		return

	def _plot_tracks(self):

		# HD LIMIT
		self.ax.plot(HD_LMT['TEFF'],HD_LMT['LOGL'], lw = HD_LMT['lw'], color = HD_LMT['c'])

		# EVOL TRACKS
		MW_tr = getTracks().MW
		for m in MW_tr:
		     if float(m) > 12 :	
 			z = np.asarray(MW_tr[m]['M_ACT']) 
 			#z = np.asarray(MW_tr[m]['TIME']/1e+6) 

    			segments = make_segments(MW_tr[m]['LTEFF'],MW_tr[m]['LOGL'])
    			lc = LineCollection(segments, array=z, cmap=self.cmap, norm=self.norm, 
						linewidth=3, alpha=1.0)
    			self.ax.add_collection(lc)
			if float(m) > 15: self.ax.text(MW_tr[m]['LTEFF'][0]+0.05,MW_tr[m]['LOGL'][0]-0.01,'%s M$_{\odot}$' % int(m), 
					size = 13, weight= 'bold')

		# LINES - THRESHOLDS
		#for l in self.llines: self.ax.axhline(y=l, color='k', ls=':')
		#for t in self.tlines: self.ax.axvline(x=np.log10(t), color='k', ls=':')

		return
