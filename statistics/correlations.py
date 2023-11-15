from plot_methods import GridTemplate
from functions import *
from constants import *
from time_series import FreqIdent
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

		self._get_data()

		return	

		self.cols_y = [c for c in self.ext_tab.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]
		super(corr_hist,self).__init__(rows_page = len(self.cols_y), cols_page = 3, row_labels = self.cols_y, 
					       fig_xlabel = PLOT_XLABEL['ls'], 
						sup_xlabels = ['Fundamental (S/N > %s)' % self.snr_show,'Fundamental', 'Harmonics'],
					       **kwargs)
		self._hist_panel()
		self.GridClose()

	def _get_data(self):

		if self.type_data == 'frequency':	

			for star in self.ext_tab['STAR']:

				freq_star = FreqIdent(star)
				star_freqs,star_snrs,star_ident=freq_star._read_freq()

				self.freqs.append(star_freqs)
				self.snrs.append(star_snrs)
				self.ident.append(star_ident)

			# fill with zero, so that column number is the same
			# enables numpy indexing/masking

			self.freqs = fill_array(self.freqs)
			self.snrs = fill_array(self.snrs)
			self.ident = fill_array(self.ident)

			#print self.freqs, self.snrs, self.ident

			return

	def _hist_panel(self):
						
		bins=np.arange(0.05,1.15,0.04);

		for prop in self.cols_y:

			ax_hsn = self.GridAx()
			ax_all = self.GridAx()
			ax_depf = self.GridAx()

			nanprop = self.ext_tab[prop].filled(np.nan)

			brp_ini = [np.nanpercentile(nanprop, 33.33),np.nanpercentile(nanprop, 66.66)] 

			tb = [-np.inf] + sorted(brp_ini) + [+np.inf]

			for i in range(len(tb)-1):

				mask_star = np.logical_and(tb[i] < self.ext_tab[prop],self.ext_tab[prop] < tb[i+1])
				mask_star = np.logical_and(mask_star,self.ext_tab[prop].mask == False)

				bin_freqs = self.freqs[mask_star]
				bin_snrs = self.snrs[mask_star]
				bin_ident = self.ident[mask_star]

				for s, ax in zip([self.snr_show,0],[ax_hsn,ax_all]):
							
					mask_fund = np.logical_and(bin_snrs > s, bin_ident == '(F)')
					hdata = bin_freqs[mask_fund]

					hdata_fit = [x for x in hdata if abs(x - np.mean(hdata)) < 3 * np.std(hdata)]
					#(mu, sigma) = norm.fit(hdata_fit)
					gparam = gamma.fit(hdata_fit, floc=0)

					if len(hdata) > 0: 
						values, _, _ =ax.hist(hdata, histtype='step', lw = 2, color = HIST_COLORS[i], 
							bins=bins, label = self._hist_legend(i,STY_LB[prop],tb[i],tb[i+1],FMT_PR[prop]))
								
						xfit = np.linspace(min(hdata), max(hdata), 100)
						area = sum(np.diff(bins) * values)

						#ax.plot(xfit,norm.pdf(xfit,mu,sigma)*area, lw=2.5, color = HIST_COLORS[i])
						ax.plot(xfit,gamma.pdf(xfit,*gparam)*area, lw=2.5, color = HIST_COLORS[i])

						ax.set_ybound(lower=0.1)
						ax.set_xbound(lower=0.01)		


				mask_depf = np.logical_and(bin_snrs > 0, bin_ident == '(H)')
				data_depf = bin_freqs[mask_depf]
				ax_depf.hist(data_depf, alpha = 0.4, histtype='step', lw = 2, color = HIST_COLORS[i],bins=bins)

			ax_hsn.text(0.7,0.3,'stbin_size %s' % mask_star.sum(), size = 8, transform=ax_hsn.transAxes)
			ax_all.legend(loc="upper right")	


	def _hist_legend(self,i,prop,low,high,frm):

		if i == 0:
			return '%s < ' % prop + frm % high
		elif i == 1:
			return frm % low + ' < %s < ' % prop + frm % high
		elif i == 2:
			return '%s > ' % prop + frm % low










