from plot_methods import GridTemplate
from functions import *
from constants import *
from time_series import FrequencyIdentifier
from astropy.table import join, Table, MaskedColumn
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
				ax.plot(jtab[c2],jtab[c1],'ko')	

				try:
				 flag_fraser = jtab['f_'+c2] == 2
				 ax.plot(jtab[c2][flag_fraser],jtab[c1][flag_fraser],'^',color='limegreen')

				 flag_haucke = jtab['f_'+c2] == 1
				 ax.plot(jtab[c2][flag_haucke],jtab[c1][flag_haucke],'ro')
				except:
				 pass
			
				ny = mask_outliers(jtab[c2], m = 4)

				if not (nx.mask.all() or ny.mask.all()):
					pears_val = pearsonr(x=nx,y=ny)[0]
					spear_val = spearmanr(x=nx,y=ny)[0]
					ax.text(0.05,0.30,'%.2f' % pears_val,size = 6,color='c',transform=ax.transAxes)
					ax.text(0.05,0.15,'%.2f' % spear_val,size = 6,color='b',transform=ax.transAxes)


			
class autocorr_scatter(corr_scatter):

	def __init__(self, vector, **kwargs):

		super(autocorr_scatter,self).__init__(ext_x = vector, ext_y = vector, auto = True, **kwargs)
		return


class tess_hist(GridTemplate):

	def __init__(self, ext_tab, snr_show = 5, **kwargs):

		self.ext_tab = ext_tab
		self.snr_show = snr_show	
		self.freqs = []
		self.snrs = []
		self.ident = []

		self.bins=np.arange(BIN_PROP['LOW'], BIN_PROP['UP'], BIN_PROP['RES'])


		len_ph = len(ext_tab)
		self.nffreq = self._masked_col(len_ph)

		self._dict_freq_tess()

		self.cols_y = [c for c in self.ext_tab.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]

		self.prop_thres = {}
		super(tess_hist,self).__init__(rows_page = len(self.cols_y), cols_page = 3, row_labels = self.cols_y, 
					       fig_xlabel = PLOT_XLABEL['ls'], 
						sup_xlabels = ['Fundamental (S/N > %s)' % self.snr_show,'Fundamental', 'Harmonics'],
					       **kwargs)
		self._hist_freq_tess()
		self.GridClose()		

		self.hist_tab = Table({'NFFREQ': self.nffreq})

		return

	def _dict_freq_tess(self):

		# Collects, classifies all TESS stored frequencies
		# Creates dictionary of stars

		self.funds = {}
		self.harms = {}
		self.combs = {}

		for s in range(len(self.ext_tab['STAR'])):

			star = self.ext_tab['STAR'][s]

			star_files = [f for f in os.listdir(TESS_LS_PATH) if 'FREQ' in f and star in f]				
			
			self.funds[star] = {}
			self.harms[star] = {}
			self.combs[star] = {}

			sectors = []
			sfunds = []; sfunds_snr = []; sfunds_ids = []
			sharms = []; sharms_snr = []; sharms_ids = []
			scombs = []; scombs_snr = []; scombs_ids = []

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
				sfunds_ids.append(sect_fid.ids[sect_fid.ff_mask & sect_fid.indiv_mask])
				sharms.append(sect_freq[sect_fid.hf_mask & sect_fid.indiv_mask])
				sharms_snr.append(sect_snr[sect_fid.hf_mask & sect_fid.indiv_mask])
				sharms_ids.append(sect_fid.ids[sect_fid.hf_mask & sect_fid.indiv_mask])
				scombs.append(sect_freq[sect_fid.cf_mask & sect_fid.indiv_mask])
				scombs_snr.append(sect_snr[sect_fid.cf_mask & sect_fid.indiv_mask])
				scombs_ids.append(sect_fid.ids[sect_fid.cf_mask & sect_fid.indiv_mask])

			self.funds[star].update({'FREQ': sfunds, 'SNR': sfunds_snr, 'IDS': sfunds_ids, 'SECT': sectors})
			self.harms[star].update({'FREQ': sharms, 'SNR': sharms_snr, 'IDS': sharms_ids, 'SECT': sectors})
			self.combs[star].update({'FREQ': scombs, 'SNR': scombs_snr, 'IDS': scombs_ids, 'SECT': sectors})

			self.nffreq[s] = len(self._convert_to_bin_averaged(np.hstack(sfunds), self.bins))
			
		freqs_to_tex_tab(d_fund = self.funds,d_harm = self.harms, d_comb = self.combs, file_='FUND_test_NEW.tex')

		return

	def _masked_col(self,len_,**mask_kwargs):
		return MaskedColumn(np.zeros(len_),mask=np.ones(len_),**mask_kwargs)

	def _hist_freq_tess(self):	

		for prop in self.cols_y:

			ax_hsn = self.GridAx()
			ax_all = self.GridAx()
			ax_harm = self.GridAx()

			nanprop = self.ext_tab[prop].filled(np.nan)
			stb_thres = [np.nanpercentile(nanprop, 33.33),np.nanpercentile(nanprop, 66.66)] 

			self.prop_thres[prop] = stb_thres

			# round TEFF in nearest 100 K
			if prop == 'TEFF' : stb_thres = [round(x,-2) for x in stb_thres]

			tb = [-np.inf] + sorted(stb_thres) + [+np.inf]

			for i in range(len(tb)-1):

				bin_mask = np.logical_and(tb[i] < self.ext_tab[prop],self.ext_tab[prop] < tb[i+1])
				bin_mask = np.logical_and(bin_mask,self.ext_tab[prop].mask == False)

				binned_stars = self.ext_tab['STAR'][bin_mask]

				binned_funds = []
				binned_funds_hsn = []
				binned_harms = []
				
				for star in binned_stars:

					# Flatten all sector data per star
					sfunds = np.hstack(self.funds[star]['FREQ'])
					sharms = np.hstack(self.harms[star]['FREQ'])
					sfunds_snr = np.hstack(self.funds[star]['SNR'])

					snr_mask = sfunds_snr > self.snr_show
					sfunds_hsn = sfunds[snr_mask]

					# Convert sector data to bin-averaged data 
					# to account max one occurence per bin per star
					binned_funds.append(self._convert_to_bin_averaged(sfunds,self.bins))
					binned_funds_hsn.append(self._convert_to_bin_averaged(sfunds_hsn,self.bins))
					binned_harms.append(self._convert_to_bin_averaged(sharms,self.bins))

				# Flatten to allow hist plotting	
				binned_funds = np.hstack(binned_funds)
				binned_funds_hsn = np.hstack(binned_funds_hsn)
				binned_harms = np.hstack(binned_harms)

				h_kwargs = {'histtype' : 'step', 'lw' : 2, 'color' : HIST_COLORS[i], 
					       'label' : self._hist_legend(i,STY_LB[prop],tb[i],tb[i+1],FMT_PR[prop]) }	

				for ax, freq_cat in zip([ax_hsn,ax_all,ax_harm],[binned_funds_hsn,binned_funds,binned_harms]):			
					self._hist_single(ax, freq_cat, self.bins, model = 'gamma', **h_kwargs)

			ax_hsn.text(0.5,0.7,'stbin_size %s' % len(binned_stars), transform=ax_hsn.transAxes)
			ax_all.legend(loc="upper right")	

		return

	def _convert_to_bin_averaged(self, data, bins):

		data_binning = np.digitize(data,bins)

		bin_averaged = []
		for k in set(data_binning):
			bin_averaged.append(np.mean([x for x,y in zip(data, data_binning) if y==k]))

		return bin_averaged
		

	def _hist_single(self, ax, data, bins, model = None, **h_kwargs):

		bd,_,_ = ax.hist(data, bins, **h_kwargs)

		if model != None :

			x_fit = np.linspace(bins[0], bins[-1], 100)
			data_fit = [x for x in data if abs(x - np.mean(data)) < 2 * np.std(data)]
			area = sum(np.diff(bins) * np.histogram(data_fit, bins)[0])

			if model == 'gamma':
				gparam = gamma.fit(data_fit, floc=0)
				ax.plot(x_fit,gamma.pdf(x_fit,*gparam)*area, lw=2.5, color = h_kwargs.get('color','k'))

			return area

		return
		
	def _hist_legend(self,i,prop,low,high,frm):

		if i == 0:
			return '%s < ' % prop + frm % high
		elif i == 1:
			return frm % low + ' < %s < ' % prop + frm % high
		elif i == 2:
			return '%s > ' % prop + frm % low
