from plot_methods import GridTemplate
import numpy as np
from functions import *
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit


class PreWhitening(GridTemplate):

	def __init__(self, lc, star, sect, save_files = False, plot_rn = True, **kwargs):

		self.lc = lc 
		self.star = star
		self.sect = sect
		self.save_files = save_files
		self.plot_rn = plot_rn
		super(PreWhitening,self).__init__(fig_xlabel='', fig_ylabel=PLOT_YLABEL['ls'], coll_x = True, params = PLOT_PARAMS['prew'], output_name= TESS_LS_PATH + '%s_%s_LSPRW' % (self.star,self.sect), col_labels = ['TESS_freq', 'TESS_time'], **kwargs)

		self.prw()	
		self.GridClose()

	def prw(self, rtol=1e-1, conv_iter=5, max_iter=50):

		lc = self.lc.remove_nans()
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

		        pow_ax = self.GridAx()
		        lc_ax = self.GridAx()
			div = make_axes_locatable(lc_ax)
			res_ax = div.append_axes("top", size="25%", pad="1%")			

			peak_power = max(pg.power)
			peak_freq = pg.frequency_at_max_power.value

			ppow_v.append(peak_power)
			pfreq_v.append(peak_freq)
			f_label_v.append('F%s' % ind)

			pow_ax.plot(pg.frequency,pg.power,'k')
			pow_ax.axvline(peak_freq,color='k',ymin=0.9,ymax=1)
			pow_ax.text(0.05,0.05,'F%s %.3f' % (ind,peak_freq),color='k', transform=pow_ax.transAxes)
			#pow_ax.hlines(max_glob,xmin=peak_freq,xmax=peak_freq+rayl,color='k')
			#pow_ax.axvline(2. / TESS_WINDOW_D,ls='-',lw=0.5)
			pow_ax.set_ylim(1.1e-6,3e-2) ; pow_ax.set_yscale('log')  
			pow_ax.set_xlim(2./TESS_WINDOW_D,25); pow_ax.set_xscale('log') 

			# Red-noise model fit
			nan_ind = np.isnan(pg.power)
			x = pg.frequency[~nan_ind] ; y = pg.power[~nan_ind]
			popt, pcov = curve_fit(self.redn,x,y)
			perr = np.sqrt(np.diag(pcov))
			SNR_pre.append(peak_power / self.redn(peak_freq,*popt))

			if self.plot_rn : pow_ax.plot(x,self.redn(x,*popt))

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

			lc_ax.plot(lc.time,lc.flux,c='k')
			lc_ax.plot(model.time,model.flux,c='r',lw=2)
			lc_ax.set_ylabel(r'$\Delta$m (mag)', size=13)
			lc_ax.invert_yaxis(); 

			prw_res = lc.copy()
			prw_res.flux = lc.flux - model.flux
			res_ax.plot(prw_res.time,prw_res.flux,'k')
			res_ax.set_ylabel('res', size=11); res_ax.yaxis.set_tick_params(labelsize=8)
			res_ax.invert_yaxis() 
			res_ax.set_xticks([])

			#lc_ax.set_ylim(0.025,-0.025); res_ax.set_ylim(0.025,-0.025)

			rn_tab["RN%s_%s" % (ind,peak_freq)] = np.append(popt,perr)
			ls_tab["P%s_%s" % (ind,peak_freq)] = pg.power 		
			model_tab["M%s_%s" % (ind,peak_freq)] = model.flux

			pg = prw_res.to_periodogram()
			model.flux = model.flux + pg.model(time=prw_res.time * u.day, frequency=pg.frequency_at_max_power).flux

			ind += 1

		ident_v = find_harmon(pfreq_v,rayl)
		ppow_v = np.array(ppow_v)
		pfreq_v = np.array(pfreq_v)
		SNR_post = ppow_v / self.redn(pfreq_v,*popt)


		if self.save_files:
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

