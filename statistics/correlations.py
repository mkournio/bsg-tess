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
		       	        cmap, norm = colorbar(self.fig, vmin=1, vmax=40, label=r'S/N',extend='max')
				for s, p, snr in zip(jtab[c1],jtab[c2], jtab['SNR_FF']):
				    flat_p = list(chain(*p))
				    flat_snr = list(chain(*snr))
				    colors = np.array([cmap(norm(k)) for k in flat_snr])
				    ax.scatter(np.array(flat_p), np.array([s] * len(flat_p)),c=colors)
			        if c2 == 'A_FF': ax.set_xlim(3E-4,5E-2)
			        ax.set_xscale('log')

                return

	def _plot_panel(self, axis, tab, k1, k2):

		e_kwargs = {'elinewidth' : 0.5, 'capsize' : 0, 'ls' : 'none'}

		mask_rot = tab['F_ROT'] == True
		mask_bright = tab['TMAG'] <= 3
		mask_contam = tab['RDMAG'] <= 2

		mask_acyg = (tab['VART'] == 'ACYG') | (tab['VART'] == 'EB+ACYG')| (tab['VART'] == 'ACYG/GCAS')
		mask_sdor = tab['VART'] == 'SDOR/L'

		axis.scatter(tab[k2][mask_rot],tab[k1][mask_rot],color='r',s=12*(7**tab['IRX'][mask_rot]))
		axis.scatter(tab[k2][~mask_rot],tab[k1][~mask_rot],color='k',s=12*(7**tab['IRX'][~mask_rot]))

		axis.plot(tab[k2][mask_bright],tab[k1][mask_bright],'cx',ms=11, mew=0.7)
		axis.plot(tab[k2][mask_contam],tab[k1][mask_contam],'co',mfc='none', mew=0.7)

		axis.plot(tab[k2][mask_acyg],tab[k1][mask_acyg],'bo', ms = 15, mfc='none')
		axis.plot(tab[k2][mask_sdor],tab[k1][mask_sdor],'gs', ms = 15, mfc='none')

		try:
		 axis.errorbar(tab[k2][mask_rot],tab[k1][mask_rot],xerr=tab['e_'+k2][mask_rot],ecolor='r',**e_kwargs)
		 axis.errorbar(tab[k2][~mask_rot],tab[k1][~mask_rot],xerr=tab['e_'+k2][~mask_rot],ecolor='k',**e_kwargs)	
		except:
		 pass	

		sig_thr = 100
		mx = mask_outliers(tab[k2], m = sig_thr)
		my = mask_outliers(tab[k1], m = sig_thr)
		if not (mx.mask.all() or my.mask.all()): spear_val = spearmanr(x=mx,y=my)[0]

		mx = mask_outliers(tab[k2][~mask_acyg], m = sig_thr)
		my = mask_outliers(tab[k1][~mask_acyg], m = sig_thr)
		if not (mx.mask.all() or my.mask.all()): spear_val_early = spearmanr(x=mx,y=my)[0]

		axis.text(0.96,0.08,r"$%.2f$" % (spear_val) + "\n" + r"$\bf{%.2f}$" % (spear_val_early),color='g',ha='right', transform=axis.transAxes)
			
		return

	def get_ax(self,row,col):

		if self.hold:
			if self.mode == 'matrix' or self.mode == 'frequency':
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
		print self.gperc

		self.outl_tab = Table()
		for c in self.x :
			self.outl_tab['OT'+c] = self._masked_col(len(data))
		self.outl_tab['STAR'] = data['STAR_1']

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

			tb = [0.999*np.nanmin(nanprop)] + sorted(stb_thres) + [1.001*np.nanmax(nanprop)]

			for col in self.y :

				ax = self.GridAx()

				i = 0
				while i < len(tb) - 1 :

					grouped_freqs = []
					grouped_acyg_freqs = []
					grouped_binary_freqs = []
					bin_mask = (tb[i] < self.data[row]) & (self.data[row] <= tb[i+1])					
					for freqs, snrs, status, star in zip(self.data[col][bin_mask],self.data['SNR_FF'][bin_mask], self.data['VART'][bin_mask], self.data['STAR_1'][bin_mask]) :
						
				    		flat_freqs = np.array(list(chain(*freqs)))
				    		flat_snrs = np.array(list(chain(*snrs)))

						mask = flat_snrs > self.snr_thres
						bin_ave_freq = convert_to_bin_averaged(flat_freqs[mask],self.bins)

						grouped_freqs.append(bin_ave_freq)
						if (status == 'ACYG') or (status == 'EB+ACYG') or (status == 'ACYG/GCAS'): 
							grouped_acyg_freqs.append(bin_ave_freq)
							#print 'ACYG: %s' % star, bin_ave_freq
						if ('EB' in status) or ('SB' in status) : 
							grouped_binary_freqs.append(bin_ave_freq)
							#print 'BINARY: %s' % star, bin_ave_freq
						#else:


					h_kwargs = {'histtype' : 'step', 'lw' : 1.5, 'color' : HIST_COLORS[i], 
				       'label' : self._hist_legend(tb[i],tb[i+1],FMT_PR[row]) }	

					grouped_freqs = list(chain(*grouped_freqs))
					model_args = self._hist_single(ax, grouped_freqs, self.bins, model = 'gamma', **h_kwargs)			
					grouped_binary_freqs = list(chain(*grouped_binary_freqs))
					if len(grouped_binary_freqs) > 0:
					 h_kwargs.update({'hatch' : 'O.', 'lw' : 0, 'label' : None})
					 self._hist_single(ax, grouped_binary_freqs, self.bins, **h_kwargs)

					grouped_acyg_freqs = list(chain(*grouped_acyg_freqs))
					if len(grouped_acyg_freqs) > 0:
				#	 h_kwargs.update({'hatch' : '//', 'alpha': 0.99, 'label' : None})
					 h_kwargs.update({ 'hatch' : None, 'lw' : 4, 'ls' : (0,(5,4)), 'label' : None})
					 self._hist_single(ax, grouped_acyg_freqs, self.bins, **h_kwargs)

					q3 = gamma(*model_args).ppf(0.75)
					q1 = gamma(*model_args).ppf(0.25) 
					iqr = q3 - q1

					for star, freqs in zip(self.data['STAR_1'][bin_mask],self.data[col][bin_mask]) :
				    		flat_freqs = np.array(list(chain(*freqs)))
						bin_ave_freq = convert_to_bin_averaged(flat_freqs,self.bins)

						#ks = kstest(bin_ave_freq, 'gamma', model_args)
						#self.pval_tab['PV'+row][self.pval_tab['STAR'] == star] = ks[1]
						self.outl_tab['OT'+row][self.outl_tab['STAR'] == star] = any(bin_ave_freq > q3 + 1.5 * iqr)
					i += 1

				lg = ax.legend(loc="upper right", fontsize=18)
				lg.set_title(STY_LB[row], prop={'size':19})	

				ax.set_yticks([5,10,15])
		return

	def _hist_single(self, ax, vals, bins, model = None, **h_kwargs):

		gparam = ()
		if model != None :

			x_fit = np.linspace(bins[0], bins[-1], 100)
			data_fit = vals # [x for x in vals if abs(x - np.mean(vals)) < 4 * np.std(vals)]
			area = sum(np.diff(bins) * np.histogram(data_fit, bins)[0])

			if model == 'gamma':
				gparam = gamma.fit(data_fit, floc=0.05)
				ax.plot(x_fit,gamma.pdf(x_fit,*gparam)*area, lw=3, color = h_kwargs.get('color','k'))

			h_kwargs['label'] += r'$\alpha$ = %.1f, $\beta$ = %.1f' % (gparam[0],1./gparam[2])

		bd,_,_ = ax.hist(vals, bins, **h_kwargs)
		ax.set_xlim(BIN_PROP['LOW']-0.03,BIN_PROP['UP'])

		return gparam
		
	def _hist_legend(self,low,high,frm):

		return r'[%s, %s]; ' % (frm % low, frm % high)


