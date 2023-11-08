from plot_methods import GridTemplate
from functions import *
import numpy as np
import lightkurve as lk
import os

class Visualize(GridTemplate):
	
	def __init__(self, level, star_ids, **kwargs):

		self.level = level
		super(Visualize,self).__init__(fig_xlabel=PLOT_XLABEL[self.level],
					       fig_ylabel=PLOT_YLABEL[self.level],**kwargs)
	
		self.rn_tab = Table(names=('STAR','w','zero','tau','gamma','e_w','e_zero','e_tau','e_gamma'), 
				    dtype=('<U16','f4', 'f4', 'f4', 'f4','f4', 'f4','f4', 'f4'))

		self._visual(star_ids)

		self.GridClose()

	def _visual(self, star_ids):

	        for star in star_ids :
			if self.level == 'lc':
				lc_files = [j for j in os.listdir(TESS_LC_PATH) if (star in j)]
				self._visual_lc(star,lc_files)

			elif self.level == 'ls':				
				rnopt = self._visual_ls(star)
				self.rn_tab.add_row(np.append(star,rnopt))
			
		return	

	def _visual_ls(self, star):

			freq_files = [j for j in os.listdir(TESS_LS_PATH) if (star in j and 'FREQS' in j) ]

			rn_chi2s = [self._get_chi2(TESS_LS_PATH+j) for j in freq_files]
			vis_sect = freq_files[np.argmin(rn_chi2s)].split('_')[1]

			ls_vis_file = '_'.join([star,vis_sect,'LSS'])
			rn_vis_file = '_'.join([star,vis_sect,'REDNOISE'])	
	
			ls_v, ls_ini, ls_end = np.loadtxt(open(TESS_LS_PATH + ls_vis_file,'rt').readlines()[:-1], delimiter=',', skiprows = 1, usecols=(0,1,-1), unpack=True)	
			rn_opt = np.loadtxt(TESS_MODEL_PATH + rn_vis_file, delimiter=',', skiprows = 1, usecols=(-1), unpack=True)	
			rnopt, rnerr = rn_opt[:4], rn_opt[4:]

			ax = self.GridAx()
			ax.plot(ls_v, ls_ini,'c'); ax.plot(ls_v, ls_end,'r')
			ax.plot(ls_v, self._redn(ls_v,*rnopt),'k--')
			ax.text(0.05,0.05,'%s (%s)' % (star,vis_sect),color='k',size=7,transform=ax.transAxes)

			ax.set_ylim(2e-6,2e-2) ; ax.set_yscale('log')  
			ax.set_xlim(3e-3,99); ax.set_xscale('log') 

			rnopt[:2] = np.log10(rnopt[:2])

			return rnopt, rnerr

	def _visual_lc(self, star, lc_files):

		if len(lc_files) > 0 :

			ax = self.GridAx()
			divider = make_axes_locatable(ax)

			ymin, ymax = self._lcs_glob_stats(lc_files)
			ym = 1.1 * max(abs(ymax),abs(ymin))

			g_lc_files = group_consec_lcs(lc_files)	
			len_gls = len(g_lc_files)

			k = 0
			for gl in g_lc_files:

				g_times = np.empty(0); g_fluxes = np.empty(0); g_sects = []
				for l in gl:
					time, flux = lc_read(TESS_LC_PATH + l)
					g_sects.append(sect_from_lc(l))	
					offset = np.nanmedian(flux)
					g_times = np.append(g_times,time)
					g_fluxes = np.append(g_fluxes,flux - offset)
				g_sect_tx = '+'.join([str(x) for x in g_sects])

				ax.plot(g_times, g_fluxes,'k'); ax.invert_yaxis(); ax.set_ylim(ym,-ym)

				ax.xaxis.set_tick_params(direction="in",labelsize=6)
				ax.yaxis.set_tick_params(labelsize=9)
				if k == 0: ax.text(0.05,0.85,star,color='r',fontsize=10,transform=ax.transAxes)
				ax.text(0.4,0.05,g_sect_tx,color='b',fontsize=9,transform=ax.transAxes	)

				if (len_gls > 1) and (k < len_gls - 1):
					d=.02
					kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
					ax.plot((1,1),(-d,+d), **kwargs) 
					ax.plot((1,1),(1-d,1+d), **kwargs) 
					ax.spines['right'].set_visible(False)

					ax = divider.append_axes("right", size="100%", pad=0.12)

					kwargs.update(transform=ax.transAxes) 
					ax.plot((0,0),(-d,+d), **kwargs) 
					ax.plot((0,0),(1-d,1+d), **kwargs) 
					ax.spines['left'].set_visible(False)
					ax.tick_params(labelleft = False) 
				k += 1

		return	

	def _get_chi2(self,path_freq):

		with open(path_freq) as f:
    			for line in f:
        			pass
    		
		return line.split()[3]

	def _redn(self,x,w,zero,tau,gamma):
		x = np.array(x)
		return w + ( zero / ( 1 + (2 * 3.14 * tau * x)**gamma))

	def _lcs_glob_stats(self,lc_files):

		minfs=[]; maxfs=[]	
		for l in lc_files:
			time, flux = np.loadtxt(TESS_LC_PATH + l, delimiter=' ', usecols=(0, 1), unpack=True)
			time_dif = time[1:] - time[:-1]; ind_nan = np.argmax(time_dif)
			flux[ind_nan-3:ind_nan+15] = np.nan; flux[:5] = np.nan; flux[-30:] = np.nan
			minfs.append(np.nanmin(flux)); maxfs.append(np.nanmax(flux))

		return min(minfs), max(maxfs)


class ProcessLevel(GridTemplate):
	
	def __init__(self, level, star_ids, folder, **kwargs):

		self.level = level
		self.LCPATH = TESS_LC_PATH+folder+'/'
		super(ProcessLevel,self).__init__(fig_xlabel=PLOT_XLABEL[self.level],
							         fig_ylabel=PLOT_YLABEL[self.level],**kwargs)		
		self._proc(star_ids)
		self.GridClose()

	def _proc(self, star_ids):
		for star in star_ids:
			lc_files = [j for j in os.listdir(self.LCPATH) if (star in j)]
			if len(lc_files) > 0: self._proc_single(star,lc_files)
		return		

	def _proc_single(self, star, lc_files):

		ax = self.GridAx()
		divider = make_axes_locatable(ax)

		if self.level == 'lc':
			ymin, ymax = self._lcs_glob_stats(lc_files)
			ym = 1.1 * max(abs(ymax),abs(ymin))

		g_lc_files = group_consec_lcs(lc_files)	
		len_gls = len(g_lc_files)

		k = 0
		for gl in g_lc_files:

			g_times = np.empty(0); g_fluxes = np.empty(0); g_sects = []
			for l in gl:
				time, flux = lc_read(self.LCPATH + l)
				g_sects.append(sect_from_lc(l))	
				offset = np.nanmedian(flux)
				g_times = np.append(g_times,time)
				g_fluxes = np.append(g_fluxes,flux - offset)
			g_sect_tx = '+'.join([str(x) for x in g_sects])

			if   self.level == 'lc':
				ax.plot(g_times, g_fluxes,'k'); ax.invert_yaxis(); ax.set_ylim(ym,-ym)
			elif self.level == 'ls':

				# Building light curve from the merged (stitched) sectors
				lc = lk.LightCurve(time=g_times,flux=g_fluxes).remove_nans()

				# Pre-whitening
				obj_prw =PreWhitening(star,g_sect_tx)
				obj_prw.prw(lc)
				obj_prw.close()

				pg = lc.to_periodogram()
				x = np.log10(np.array(pg.frequency))
				y = np.log10(np.array(pg.power))

				ax.plot(x,y,'b')
				ax.set_ylim(-5.9,-2.1)	


			ax.xaxis.set_tick_params(direction="in",labelsize=6)
			ax.yaxis.set_tick_params(labelsize=9)
			if k == 0: ax.text(0.05,0.85,star,color='r',fontsize=10,transform=ax.transAxes)
			ax.text(0.4,0.05,g_sect_tx,color='b',fontsize=9,transform=ax.transAxes	)

			if (len_gls > 1) and (k < len_gls - 1):
				d=.02
				kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
				ax.plot((1,1),(-d,+d), **kwargs) 
				ax.plot((1,1),(1-d,1+d), **kwargs) 
				ax.spines['right'].set_visible(False)

				ax = divider.append_axes("right", size="100%", pad=0.12)

				kwargs.update(transform=ax.transAxes) 
				ax.plot((0,0),(-d,+d), **kwargs) 
				ax.plot((0,0),(1-d,1+d), **kwargs) 
				ax.spines['left'].set_visible(False)
				ax.tick_params(labelleft = False) 
					

			k += 1

		return	

	def _redn(self,x,w,zero,tau,gamma):
		x = np.array(x)
		return w + ( zero / ( 1 + (2 * 3.14 * tau * x)**gamma))

	def _lcs_glob_stats(self,lc_files):

		minfs=[]; maxfs=[]	
		for l in lc_files:
			time, flux = np.loadtxt(self.LCPATH + l, delimiter=' ', usecols=(0, 1), unpack=True)
			time_dif = time[1:] - time[:-1]; ind_nan = np.argmax(time_dif)
			flux[ind_nan-3:ind_nan+15] = np.nan; flux[:5] = np.nan; flux[-30:] = np.nan
			minfs.append(np.nanmin(flux)); maxfs.append(np.nanmax(flux))

		return min(minfs), max(maxfs)

