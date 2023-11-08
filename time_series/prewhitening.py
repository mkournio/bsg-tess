from plot_methods import GridTemplate
import numpy as np
from functions import *
import pandas as pd
import lightkurve as lk
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import os 
import matplotlib.ticker
from scipy.optimize import curve_fit



class PreWhitening(object):

	def __init__(self, star, sect, store_data = True):

		self.star = star
		self.sect = sect
		self.store_data = store_data
		self.pdf = PdfPages(TESS_LS_PATH + '%s_%s_LSPRW.pdf' % (self.star,self.sect))

		return

	def prw(self, lc, rtol=1e-1, conv_iter=5, max_iter=50):

		lc = lc.remove_nans()
		prw_lc = lc.copy()
		pg = prw_lc.to_periodogram()
		model = pg.model(time=prw_lc.time * u.day,frequency=pg.frequency_at_max_power)
		max_glob = max(pg.power)

		# Initializing tabs
		model_tab = pd.DataFrame({"DATE":lc.time, "FLUX":lc.flux})
		ls_tab = pd.DataFrame({"FREQUENCY":pg.frequency})
		rn_tab = pd.DataFrame()
		pfreq_v = []; ppow_v = [];  f_label_v=[]
		CHI2_vect = []
		SNR_pre = []

		# Initializing variables/constants
		rayl = 1 / (lc.time[-1] - lc.time[0])
		term_index = 0
		ind = 0

		while (term_index < conv_iter) and (ind < max_iter) : 		

			if ind % LSPRW_PLOT_NUM == 0:
				self.fig = plt.figure(figsize=(10,12))
				gs = GridSpec(LSPRW_PLOT_NUM, 2, figure=self.fig, hspace=0.3, wspace=0.4)

			plot_ind = ind % LSPRW_PLOT_NUM
			ax1 = self.fig.add_subplot(gs[plot_ind, 0])
			ax2 = self.fig.add_subplot(gs[plot_ind, 1])
			div = make_axes_locatable(ax2)
			ax3 = div.append_axes("bottom", "40%", pad=0.)

			peak_power = max(pg.power)
			peak_freq = pg.frequency_at_max_power.value

			ppow_v.append(peak_power)
			pfreq_v.append(peak_freq)
			f_label_v.append('F%s' % ind)

			pg.plot(ax=ax1); ax1.axvline(peak_freq,color='r',ymin=0.9,ymax=1); ax1.set_ylim(1e-6,None) ; ax1.set_yscale('log')  
			ax1.set_xlim(3e-3,100); ax1.set_xscale('log') 
			ax1.text(0.05,0.05,'F%s %s' % (ind,peak_freq),color='r',size=7,transform=ax1.transAxes)
			ax1.hlines(max_glob,xmin=peak_freq,xmax=peak_freq+rayl,color='r')
			ax1.axvline(2*rayl,ls='-',lw=0.3)

			# Red-noise model fit
			nan_ind = np.isnan(pg.power)
			x = pg.frequency[~nan_ind] ; y = pg.power[~nan_ind]
			popt, pcov = curve_fit(self.redn,x,y)
			perr = np.sqrt(np.diag(pcov))
			SNR_pre.append(peak_power / self.redn(peak_freq,*popt))

			ax1.plot(x,self.redn(x,*popt))

			CHI2 =  np.sum( ((self.redn(x,*popt)-y)**2 ) / self.redn(x,*popt))
			CHI2_vect.append(CHI2);
			try:
				rel_change = abs(CHI2_old - CHI2) / CHI2_old
				if rel_change < rtol:
					term_index += 1
				else:
					term_index = 0
			except:
				pass
			CHI2_old = CHI2

			lc.plot(ax=ax2,c='k'); model.plot(ax=ax2,c='r',ylabel=r'$\Delta$m (mag)',label='')
			ax2.invert_yaxis(); 

			prw_res = lc.copy()
			prw_res.flux = lc.flux - model.flux
			prw_res.plot(ax=ax3,c='k', ylabel='residuals'); ax3.invert_yaxis() 

			rn_tab["RN%s_%s" % (ind,peak_freq)] = np.append(popt,perr)
			ls_tab["P%s_%s" % (ind,peak_freq)] = pg.power 		
			model_tab["M%s_%s" % (ind,peak_freq)] = model.flux

			pg = prw_res.to_periodogram()
			model.flux = model.flux + pg.model(time=prw_res.time * u.day, frequency=pg.frequency_at_max_power).flux

			ind += 1
			if ind % LSPRW_PLOT_NUM == 0: 
				self.pdf.savefig(self.fig)
				self.fig.clf(); plt.close(self.fig)


		if ind % LSPRW_PLOT_NUM != 0: self.pdf.savefig(self.fig)

		ident_v = find_harmon(pfreq_v,rayl)

		ppow_v = np.array(ppow_v)
		pfreq_v = np.array(pfreq_v)
		SNR_post = ppow_v / self.redn(pfreq_v,*popt)


		if self.store_data:
			ls_tab.to_csv(TESS_LS_PATH + self._head_format('LSS'), index=False)
			model_tab.to_csv(TESS_MODEL_PATH + self._head_format('SYNTHLC'), index=False)
			rn_tab.to_csv(TESS_MODEL_PATH + self._head_format('REDNOISE'), index=False)

			with open(TESS_LS_PATH + self._head_format('FREQS'),'w') as freq_f:
				freq_f.write('Rayleigh %s\n' % rayl)
				for i,j,k,c,l,m,n in zip(f_label_v,pfreq_v,ppow_v,CHI2_vect,SNR_pre,SNR_post,ident_v):
    					freq_f.write("%s %.5f %.3E %.3f %.2f %.2f %s\n" % (i,j,k,c,l,m,n))
			freq_f.close()

		return 

	def _head_format(self,id):
		return '%s_%s_%s' % (self.star,self.sect,id)

	def redn(self,x,w,zero,tau,gamma):
		x = np.array(x)
		return w + ( zero / ( 1 + (2 * 3.14 * tau * x)**gamma))

	def close(self):

		self.fig.clf(); plt.close(self.fig)
		self.pdf.close()

		return

