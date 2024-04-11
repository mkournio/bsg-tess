from plot_methods import *
from functions import *
from constants import *
from time_series import FrequencyIdentifier
from astropy.table import join, Table, MaskedColumn
import numpy as np
from scipy.stats import norm, gamma, kstest
from scipy.stats.mstats import pearsonr, spearmanr, describe
from itertools import chain
import pandas as pd

class corr_scatter(GridTemplate):

	def __init__(self, data, x, y, mode = 'matrix', hold = False, **kwargs):

		jtab = data.copy()
		self.x = x
		self.y = y
		self.hold = hold
		self.mode = mode
		self._validate()

		mode_grid = {}
		if self.mode != 'zip':
			mode_grid = {'rows_page' : len(self.x), 'cols_page' : len(self.y), 'coll_x' : True, 'coll_y' : True,
				     'row_labels' : self.x, 'col_labels' : self.y} 
		super(corr_scatter,self).__init__(fig_xlabel='', fig_ylabel='', params = PLOT_PARAMS['cr'], **dict(kwargs,**mode_grid))

		self._plotting(jtab)

		if not hold : 
			self.GridClose()

	def _validate(self):

		if self.mode == 'zip' and len(self.x) != len(self.y):
			raise IndexError('Zip mode is activated but key vectors have not same size.')
		pass

	def _plotting(self, jtab):

		if self.mode == 'zip':

			for c1, c2 in zip(self.x, self.y) :
				ax = self.GridAx()
				self._plot_panel(ax,jtab,c1,c2)

				ax.set_xlabel(STY_LB[c2])
				ax.set_ylabel(STY_LB[c1])

		elif self.mode == 'matrix':

			for c1 in self.x :

			 prl_y = np.nanmax(jtab[c1]) - np.nanmin(jtab[c1])
			 prm_y = 0.5*(np.nanmax(jtab[c1])+np.nanmin(jtab[c1]))
			 min_yaxis = prm_y - 0.6*prl_y 
			 max_yaxis = prm_y + 0.6*prl_y

			 for c2 in self.y :

			 	prl_x = np.nanmax(jtab[c2]) - np.nanmin(jtab[c2])
			 	prm_x = 0.5*(np.nanmax(jtab[c2])+np.nanmin(jtab[c2]))
			 	min_xaxis = prm_x - 0.6*prl_x 
			 	max_xaxis = prm_x + 0.6*prl_x

				ax = self.GridAx()
				self._plot_panel(ax,jtab,c1,c2)

				ax.set_xlim(min_xaxis,max_xaxis)
				ax.set_ylim(min_yaxis,max_yaxis)

		elif self.mode == 'frequency':

			for c1 in self.x :
			 for c2 in self.y :

				ax = self.GridAx()
		       	        cmap, norm = colorbar(self.fig, vmin=1, vmax=50, label=r'S/N',extend='max')
				for s, p, snr in zip(jtab[c1],jtab[c2], jtab['SNR_FF']):

				    flat_p = list(chain(*p))
				    flat_snr = list(chain(*snr))
				    colors = np.array([cmap(norm(k)) for k in flat_snr])
				    ax.scatter(np.array(flat_p), np.array([s] * len(flat_p)),c=colors)
                                ax.set_xscale('log')		
				

                return

	def _plot_panel(self, axis, tab, k1, k2):

		e_kwargs = {'elinewidth' : 0.5, 'capsize' : 0, 'ls' : 'none'}

		try:
		 flag_fraser = tab['f_'+k1] == 2
		 axis.plot(tab[k2][flag_fraser],tab[k1][flag_fraser],'g^')
		 axis.errorbar(tab[k2][flag_fraser],tab[k1][flag_fraser],xerr=tab['e_'+k2][flag_fraser],ecolor='g',**e_kwargs)
	
		 flag_haucke = tab['f_'+k1] == 1
		 axis.plot(tab[k2][flag_haucke],tab[k1][flag_haucke],'bs')
		 axis.errorbar(tab[k2][flag_haucke],tab[k1][flag_haucke],xerr=tab['e_'+k2][flag_haucke],ecolor='r',**e_kwargs)	
		except:
		 axis.plot(tab[k2],tab[k1],'ko')
		 try:
		  axis.errorbar(tab[k2],tab[k1],xerr=tab['e_'+k2],ecolor='k',**e_kwargs)	
		 except:
		  pass
	
		mx = mask_outliers(tab[k1], m = 4)
		my = mask_outliers(tab[k2], m = 4)

		if not (mx.mask.all() or my.mask.all()):
			spear_val = spearmanr(x=mx,y=my)[0]
			axis.text(0.70,0.05,'%.2f' % spear_val,color='r',transform=axis.transAxes)
			
		return

	def get_ax(self,row,col):

		if self.hold:
			if self.mode == 'matrix':
				return self.fig.axes[row*len(self.y)+col]
			elif self.mode == 'zip':
				return self.fig.axes[row*self.cols_page+col]
		else:
			print 'Warning: Figure has closed. Use hold = True'
			pass
		return


class freq_hist(GridTemplate):

	def __init__(self, data, x, y, snr_thres = 5, ngroups = 3, **kwargs):

		self.data = data
		self.x = x
		self.y = y
		self.snr_thres = snr_thres	
		self.freqs = []
		self.snrs = []

		self.gperc = np.linspace(0,100,ngroups+1)[1:-1]

		self.pval_tab = Table()
		for c in self.x :
			self.pval_tab['PV'+c] = self._masked_col(len(data))
			self.pval_tab['OT'+c] = self._masked_col(len(data))
		self.pval_tab['STAR'] = data['STAR_1']

		self.bins=np.arange(BIN_PROP['LOW'], BIN_PROP['UP'], BIN_PROP['RES'])


		self.thres_group = {}

		super(freq_hist,self).__init__(rows_page = len(self.x), cols_page = len(self.y), col_labels = self.y, coll_x = True, fig_xlabel='', fig_ylabel='#stars', params = PLOT_PARAMS['ls'], **kwargs)

		self._hist_freq()

		self.GridClose()		

		return


	def _masked_col(self,len_,**mask_kwargs):
		return MaskedColumn(np.zeros(len_),mask=np.ones(len_),**mask_kwargs)

	def _hist_freq(self):	

		for row in self.x:

			nanprop = self.data[row].filled(np.nan)
			
			stb_thres = [np.nanpercentile(nanprop, x) for x in self.gperc]
			self.thres_group[row] = stb_thres

			tb = [0.9*np.nanmin(nanprop)] + sorted(stb_thres) + [1.1*np.nanmax(nanprop)]

			for col in self.y :

				ax = self.GridAx()

				i = 0
				while i < len(tb) - 1 :

					grouped_freqs = []
					bin_mask = (tb[i] < self.data[row]) & (self.data[row] <= tb[i+1])

					for freqs, snrs in zip(self.data[col][bin_mask],self.data['SNR_FF'][bin_mask]) :
						
				    		flat_freqs = np.array(list(chain(*freqs)))
				    		flat_snrs = np.array(list(chain(*snrs)))

						mask = flat_snrs > self.snr_thres
						bin_ave_freq = self._convert_to_bin_averaged(flat_freqs[mask],self.bins)
						grouped_freqs.append(bin_ave_freq)

					h_kwargs = {'histtype' : 'step', 'lw' : 2, 'color' : HIST_COLORS[i], 
				       'label' : self._hist_legend(tb[i],tb[i+1],FMT_PR[row]) }	
					model_args = self._hist_single(ax, list(chain(*grouped_freqs)), self.bins, model = 'gamma', **h_kwargs)
					q3 = gamma(*model_args).ppf(0.75)
					q1 = gamma(*model_args).ppf(0.25) 
					iqr = q3 - q1

					for star, freqs in zip(self.data['STAR_1'][bin_mask],self.data[col][bin_mask]) :
				    		flat_freqs = np.array(list(chain(*freqs)))
						bin_ave_freq = self._convert_to_bin_averaged(flat_freqs,self.bins)

						ks = kstest(bin_ave_freq, 'gamma', model_args)
						self.pval_tab['PV'+row][self.pval_tab['STAR'] == star] = ks[1]
						self.pval_tab['OT'+row][self.pval_tab['STAR'] == star] = any(bin_ave_freq > q3 + 1.5 * iqr)

					i += 1

				ax.legend(title=STY_LB[row], loc="upper right")	
		return

	def _convert_to_bin_averaged(self, data, bins):

		data_binning = np.digitize(data,bins)

		bin_averaged = []
		for k in set(data_binning):
			bin_averaged.append(np.mean([x for x,y in zip(data, data_binning) if y==k]))

		return bin_averaged
		

	def _hist_single(self, ax, vals, bins, model = None, **h_kwargs):

		if model != None :

			x_fit = np.linspace(bins[0], bins[-1], 100)
			data_fit = [x for x in vals if abs(x - np.mean(vals)) < 3 * np.std(vals)]
			area = sum(np.diff(bins) * np.histogram(data_fit, bins)[0])

			if model == 'gamma':
				gparam = gamma.fit(data_fit, floc=0.06)
				ax.plot(x_fit,gamma.pdf(x_fit,*gparam)*area, lw=2.5, color = h_kwargs.get('color','k'))

			h_kwargs['label'] += r'$\alpha$ = %.1f, $\beta$ = %.1f' % (gparam[0],1./gparam[2])

		bd,_,_ = ax.hist(vals, bins, **h_kwargs)
		ax.set_xlim(BIN_PROP['LOW']-0.03,BIN_PROP['UP']-0.03)

		return gparam
		
	def _hist_legend(self,low,high,frm):

		return r'[%s, %s]; ' % (frm % low, frm % high)


