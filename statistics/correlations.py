from plot_methods import GridTemplate
from functions import *
from constants import *
from time_series import FrequencyIdentifier
from astropy.table import join
import numpy as np
from scipy.stats import norm, gamma 
from scipy.stats.mstats import pearsonr, spearmanr, describe

class corr_scatter(GridTemplate):

	def __init__(self, ext_x, ext_y, match_keys = 'STAR', auto = False, **kwargs):

		self.ext_x = ext_x
		self.ext_y = ext_y
		self.match_keys = match_keys
		self.auto = auto

		self.cols_x = [c for c in self.ext_x.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]
		self.cols_y = [c for c in self.ext_y.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]

		super(corr_scatter,self).__init__(rows_page = len(self.cols_x), cols_page = len(self.cols_y),
						 row_labels = self.cols_x, col_labels = self.cols_y, **kwargs)

		self._scatter_panel()
		self.GridClose()	

	def _scatter_panel(self):

		if self.auto:
			jtab = self.ext_x
		else:
			jtab = join(self.ext_x, self.ext_y, keys=self.match_keys)

		if 'TEFF' in jtab.columns: jtab['TEFF'] = np.log10(jtab['TEFF'])

		for c1 in self.cols_x :

			nx = mask_outliers(jtab[c1], m = 4)

			for c2 in self.cols_y :

				ax = self.GridAx()

				flag_fraser = jtab['f_'+c2] == 2
				ax.plot(jtab[c2][flag_fraser],jtab[c1][flag_fraser],'^',color='limegreen')

				flag_haucke = jtab['f_'+c2] == 1
				ax.plot(jtab[c2][flag_haucke],jtab[c1][flag_haucke],'ro')				

				ny = mask_outliers(jtab[c2], m = 4)

				ax.plot(ny,nx,'k.',markersize=2)

				if not (nx.mask.all() or ny.mask.all()):
					pears_val = pearsonr(x=nx,y=ny)[0]
					spear_val = spearmanr(x=nx,y=ny)[0]
					ax.text(0.05,0.30,'%.2f' % pears_val,size = 6,color='c',transform=ax.transAxes)
					ax.text(0.05,0.15,'%.2f' % spear_val,size = 6,color='b',transform=ax.transAxes)

				#ax.plot(jtab[c2],jtab[c1],'k.')	
			
class autocorr_scatter(corr_scatter):

	def __init__(self, vector, **kwargs):

		super(autocorr_scatter,self).__init__(ext_x = vector, ext_y = vector, auto = True, **kwargs)
		return


class corr_hist(GridTemplate):

	def __init__(self, ext_tab, snr_show = 5, type_data = 'frequency', **kwargs):

		self.ext_tab = ext_tab
		self.type_data = type_data
		self.snr_show = snr_show	
		self.freqs = []
		self.snrs = []
		self.ident = []

		self.bins=np.arange(BIN_PROP['LOW'], BIN_PROP['UP'], BIN_PROP['RES'])

		self._dict_freq_tess()

		self.cols_y = [c for c in self.ext_tab.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]
		super(corr_hist,self).__init__(rows_page = len(self.cols_y), cols_page = 3, row_labels = self.cols_y, 
					       fig_xlabel = PLOT_XLABEL['ls'], 
						sup_xlabels = ['Fundamental (S/N > %s)' % self.snr_show,'Fundamental', 'Harmonics'],
					       **kwargs)
		self._hist_freq_tess()
		self.GridClose()

	def _dict_freq_tess(self):

		self.funds = {}
		self.harms = {}
		self.combs = {}

		for star in self.ext_tab['STAR']:

			star_files = [f for f in os.listdir(TESS_LS_PATH) if 'FREQ' in f and star in f]				
			
			self.funds[star] = {}
			self.harms[star] = {}
			self.combs[star] = {}

			sectors = []
			sfunds = []; sfunds_snr = []
			sharms = []; sharms_snr = []
			scombs = []; scombs_snr = []

			for sf in star_files:

				sect_freq = []
				sect_snr = []

				sectors.append(sf.split('_')[1])

				with open(TESS_LS_PATH+sf) as file_:
					first_line = file_.readline()
					rayleigh = float(first_line.split()[1])
					for line in file_:
				     		if float(line.split()[1]) > 2 * rayleigh:
							fr, snr = map(line.split().__getitem__,[1,5])
							sect_freq.append(float(fr))
							sect_snr.append(float(snr))					
				file_.close()

				sect_freq = np.array(sect_freq)
				sect_snr = np.array(sect_snr)

				m_ind = np.argmin(sect_freq[:2])
				sect_freq[0], sect_freq[m_ind] = sect_freq[m_ind], sect_freq[0]
				sect_snr[0], sect_snr[m_ind] = sect_snr[m_ind], sect_snr[0]

				sect_fid = FrequencyIdentifier(sect_freq, resolution = rayleigh)

				sfunds.append(sect_freq[sect_fid.ff_mask & sect_fid.indiv_mask])
				sfunds_snr.append(sect_snr[sect_fid.ff_mask & sect_fid.indiv_mask])
				sharms.append(sect_freq[sect_fid.hf_mask & sect_fid.indiv_mask])
				sharms_snr.append(sect_snr[sect_fid.hf_mask & sect_fid.indiv_mask])
				scombs.append(sect_freq[sect_fid.cf_mask & sect_fid.indiv_mask])
				scombs_snr.append(sect_snr[sect_fid.cf_mask & sect_fid.indiv_mask])

			self.funds[star].update({'FREQ': sfunds, 'SNR': sfunds_snr, 'SECT': sectors})
			self.harms[star].update({'FREQ': sharms, 'SNR': sharms_snr, 'SECT': sectors})
			self.combs[star].update({'FREQ': scombs, 'SNR': scombs_snr, 'SECT': sectors})

		return

	def _hist_freq_tess(self):	

		for prop in self.cols_y:

			ax_hsn = self.GridAx()
			ax_all = self.GridAx()
			ax_harm = self.GridAx()

			nanprop = self.ext_tab[prop].filled(np.nan)
			brp_ini = [np.nanpercentile(nanprop, 33.33),np.nanpercentile(nanprop, 66.66)] 

			tb = [-np.inf] + sorted(brp_ini) + [+np.inf]

			for i in range(len(tb)-1):

				bin_mask = np.logical_and(tb[i] < self.ext_tab[prop],self.ext_tab[prop] < tb[i+1])
				bin_mask = np.logical_and(bin_mask,self.ext_tab[prop].mask == False)

				binned_stars = self.ext_tab['STAR'][bin_mask]

				binned_funds = []
				binned_funds_hsn = []
				binned_harms = []
				
				for star in binned_stars:

					sfunds = np.hstack(self.funds[star]['FREQ'])
					sharms = np.hstack(self.harms[star]['FREQ'])

					sfunds_snr = np.hstack(self.funds[star]['SNR'])
					snr_mask = sfunds_snr > self.snr_show

					sfunds_hsn = sfunds[snr_mask]
				
					for k in set(np.digitize(sfunds,self.bins)):
						mask = np.digitize(sfunds,self.bins) == k
						binned_funds.append(np.mean([x for x,y in zip(sfunds,mask) if y]))
					for k in set(np.digitize(sfunds_hsn,self.bins)):
						mask = np.digitize(sfunds_hsn,self.bins) == k
						binned_funds_hsn.append(np.mean([x for x,y in zip(sfunds_hsn,mask) if y]))
					for k in set(np.digitize(sharms,self.bins)):
						mask = np.digitize(sharms,self.bins) == k
						binned_harms.append(np.mean([x for x,y in zip(sharms,mask) if y]))
			
				h_kwargs = {'histtype' : 'step', 'lw' : 2, 'color' : HIST_COLORS[i], 
					       'label' : self._hist_legend(i,STY_LB[prop],tb[i],tb[i+1],FMT_PR[prop]) }
				

				self._hist_single(ax_hsn, binned_funds_hsn, self.bins, model = 'gamma', **h_kwargs)
				self._hist_single(ax_all, binned_funds, self.bins, model = 'gamma', **h_kwargs)
				self._hist_single(ax_harm, binned_harms, self.bins, model = 'gamma', **h_kwargs)

			#ax_hsn.text(0.7,0.3,'stbin_size %s' % mask_star.sum(), size = 8, transform=ax_hsn.transAxes)
			ax_all.legend(loc="upper right")	

		return

	def _hist_single(self, ax, data, bins, model = None, **h_kwargs):

		bd,_,_ = ax.hist(data, bins, **h_kwargs)

		if model != None :

			area = sum(np.diff(bins) * bd)
			x_fit = np.linspace(bins[0], bins[-1], 100)
			data_fit = [x for x in data if abs(x - np.mean(data)) < 2 * np.std(data)]

			if model == 'gamma':
				gparam = gamma.fit(data_fit, floc=0)
				ax.plot(x_fit,gamma.pdf(x_fit,*gparam)*area, lw=2.5, color = h_kwargs.get('color','k'))

		return
		
	def _hist_legend(self,i,prop,low,high,frm):

		if i == 0:
			return '%s < ' % prop + frm % high
		elif i == 1:
			return frm % low + ' < %s < ' % prop + frm % high
		elif i == 2:
			return '%s > ' % prop + frm % low
