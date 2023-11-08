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

class ExtLCs(GridTemplate):

	def __init__(self,stars_id,stars_ra,stars_dec,stars_tic,plots_per_grid,**kwargs):

		fl_id = 0 
		while os.path.exists(TESS_LC_PATH+'EXTLC_%s' % fl_id): fl_id += 1
		self.EXTLCpath = TESS_LC_PATH+'EXTLC_%s' % fl_id
		os.mkdir(self.EXTLCpath)

		self.fSPOC = os.listdir(TESS_SPOC)
		self.fFFI = os.listdir(TESS_TPF_PATH)
		super(ExtLCs,self).__init__(plots_per_grid, cols_page = PLOT_XLC_NCOL, fig_xlabel=PLOT_LC_XLABEL,**kwargs)	
		self._extract_LCs(stars_id, stars_ra, stars_dec,stars_tic,**kwargs)
		self.GridClose()

	def _extract_LCs(self, stars_id, stars_ra, stars_dec, stars_tic, **kwargs) :

		    for s in range(len(stars_id)) :

			star = stars_id[s]
			ra = stars_ra[s]
			dec = stars_dec[s]
			tic = stars_tic[s]

			if 'thmask' in kwargs:
				thmask = kwargs["thmask"][s]
			else:
				thmask = FFI_MASK			

			spoc_files = [f for f in self.fSPOC if str(tic) in f]
			spoc_sects = [int(f.split('-')[1][1:]) for f in spoc_files]

			ffi_files = [j for j in self.fFFI if (str(ra)[:5] in j) and (str(dec)[:5] in j)]
			ffi_sects = [int(f.split('-')[1][1:]) for f in ffi_files]

			try:
				sspoc = kwargs["sspoc"][s]; sffi = kwargs["sffi"][s]
				sspoc = sect_unmasked(sspoc); sffi = sect_unmasked(sffi);
				spoc_files, spoc_sects = [x for x, y in zip(spoc_files,spoc_sects) if y in sspoc], sspoc
				ffi_files, ffi_sects = [x for x, y in zip(ffi_files,ffi_sects) if y in sffi], sffi
			except:
				pass

			star_sects =  sorted(set(spoc_sects + ffi_sects))
			if len(star_sects) > 0: self._extract_single(star,ra,dec,star_sects,spoc_files,ffi_files,thmask)

		    return

 	def _extract_single(self, star, ra, dec, star_sects, spoc_files, ffi_files, thmask):

		sect_old = -2 
		for sect in star_sects:

			axes = self.GridAx()

			spoc_sect = [f for f in spoc_files if int(f.split('-')[1][1:]) == sect]
			ffi_sect = [f for f in ffi_files if int(f.split('-')[1][1:]) == sect]
			pca=1; ndeg=2

			if bool(spoc_sect):
				print('%s, sector %s: SPOC lightcurve processing..' % (star,sect))
				spoc_lc = lk.TessLightCurveFile(TESS_SPOC + spoc_sect[0]).get_lightcurve('PDCSAP_FLUX')

				x_act = spoc_lc.time
				y_act = spoc_lc.flux
				yerr_act = spoc_lc.flux_err
				ltype = 'SPOC' 		

			elif bool(ffi_sect):
			  try:
				print('%s, sector %s: FFI, extracting and processing..' % (star,sect))
				tpf = lk.TessTargetPixelFile(TESS_TPF_PATH + ffi_sect[0])
				bk_mask = ~tpf.create_threshold_mask(THR_BCKG, reference_pixel=None)
				bk_lc = tpf.to_lightcurve(aperture_mask=bk_mask).remove_nans()
				RA_pix, DE_pix, gaia_sizes, gaia_ind, star_row, star_col  = self._gaia(tpf,ra,dec)

				mask = getmask(tpf,star_row,star_col,thres=thmask)
				save_mask(mask,TESS_MASK_PATH+'%s_%s_MASK' % (star,tpf.sector))

				tpf.plot(ax=axes[0],aperture_mask=mask,mask_color='#FD110D',show_colorbar=False,\
					bkg_mask=bk_mask,bkg_mask_color='w',title='')
				axes[0].scatter(RA_pix,DE_pix,s=gaia_sizes, marker='.', c='c') 
				axes[0].text(0.05,0.90,'%s' % thmask,color='c',size=12,transform=axes[0].transAxes)

				if gaia_ind != None: axes[0].scatter(RA_pix[gaia_ind],DE_pix[gaia_ind],s=GAIA_UPMARK,\
							marker='x',color='c',linewidths=2)

				raw_lc = tpf.to_lightcurve(aperture_mask=mask).remove_nans()
				corrected_lc1 = lccor(tpf, mask, bk_mask, pca)
				x_act = corrected_lc1.time
				y_act = corrected_lc1.flux
				yerr_act = corrected_lc1.flux_err
				ltype = 'FFI'
				
		          except:
				print 'ERROR IN FFI'
			   	pass	

			if sect == sect_old + 1:
				offset = ymed_old - np.nanmedian(y_act)
				y_act += offset

			x_n, yfit,dm,e_dm= fit_pol(x_act,y_act,yerr_act,ndeg,unit='mag')

			axes[1].plot(x_act,y_act,LC_COLOR[ltype])
			axes[1].plot(x_n,yfit,LC_COLOR['polyfit'])
			axes[2].plot(x_n,dm,LC_COLOR[ltype])

			save_three_col(x_n,dm,e_dm,self.EXTLCpath+'/%s_%s_DMAG_%s' % (star,sect,ltype))

			sect_old = sect
			ymed_old = np.nanmedian(y_act)

			axes[0].set_ylabel('%s (%s)' % (star,sect), fontsize=16); axes[0].set_xlabel('')
			axes[2].invert_yaxis()			


		return

	def _gaia(self,tpf,ra,dec):

			Gaia3_tpf = Gaia(tpf,cat='Gaia3')
			RA_pix, DE_pix, BPRP, RP, Ga_IDs = Gaia3_tpf.get_prop()
			gaia_ind = None

			try:
				from astroquery.vizier import Vizier
				from astropy.coordinates import SkyCoord
				from astropy import units as u

				GA_query = Vizier(catalog='I/355/gaiadr3').query_region(SkyCoord(ra=ra,dec=dec,unit=(u.deg,u.deg),\
					 frame='icrs'), radius = 10 * u.arcsec)[0]; GA_star =  GA_query['Source'][0]
				gaia_ind = np.where(Ga_IDs == GA_star)
				star_row = int(DE_pix[gaia_ind] - tpf.row); star_col = int(RA_pix[gaia_ind] - tpf.column)
			except:
				star_row = int(0.5*tpf.shape[2]); star_col = int(0.5*tpf.shape[1])

			rp_min = np.where(RP == min(RP))
			RPdiff =  RP - RP[rp_min]
			gaia_sizes = GAIA_UPMARK / 2**RPdiff

			return RA_pix, DE_pix, gaia_sizes, gaia_ind, star_row, star_col

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

class Correlation(GridTemplate):

	def __init__(self, ext_x, ext_y, match_keys, **kwargs):

		self.ext_x = ext_x
		self.ext_y = ext_y
		self.match_keys = match_keys

		self.cols_x = [c for c in self.ext_x.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]
		self.cols_y = [c for c in self.ext_y.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]

		super(Correlation,self).__init__(rows_page = len(self.cols_x), cols_page = len(self.cols_y),
						 row_labels = self.cols_x, col_labels = self.cols_y, **kwargs)

		self._corr_panel()
		self.GridClose()	

	def _corr_panel(self):

		from astropy.table import join
		jtab = join(self.ext_x, self.ext_y, keys=self.match_keys)

		try:
			jtab['TEFF'] = np.log10(jtab['TEFF'])
		except:
			pass

		for c1 in self.cols_x :
			for c2 in self.cols_y :

				ax = self.GridAx()

				flag_fraser = jtab['f_'+c2] == 3
				ax.plot(jtab[c2][flag_fraser],jtab[c1][flag_fraser],'g^')

				flag_haucke = jtab['f_'+c2] == 1
				ax.plot(jtab[c2][flag_haucke],jtab[c1][flag_haucke],'ro')

				#ax.plot(jtab[c2],jtab[c1],'k.')				




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



#if __name__ == '__main__':

import argparse
from astropy.table import Table
from astroquery.xmatch import XMatch
from astropy.io import ascii

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename",default= None, help="File containing the star table")
parser.add_argument("-l1","--from_line", type=int, default = 1, help="First reading line")
parser.add_argument("-l2","--to_line", type=int, default = None, help="End reading line")	
args = parser.parse_args()
	
	#### READ TABLE
inputf = ascii.read(args.filename, header_start=0, delimiter=",",data_start=args.from_line, data_end=args.to_line)
inputf['RA'], inputf['DEC'] = to_deg(inputf['RA'],inputf['DEC'])

coord = Table({'STAR' : inputf['STAR'], 'RA' : inputf['RA'], 'DEC' : inputf['DEC'], 
		       'SSPOC': inputf['SSPOC'], 'SFFI' : inputf['SFFI'], 'MASK' : inputf['MASK']})


data = XMatch.query(cat1=coord, cat2='vizier:IV/39/tic82',  max_distance=2*u.arcsec, colRA1='RA', colDec1='DEC')
fdata =  match_tabs(coord,data)
fdata.pprint(max_lines=-1)
	
	#### EXTRACT LIGHTCURVES (ALL)
	#ExtLCs(data['STAR'],data['RA'],data['DEC'],data['TIC'],4,pdf_name=args.filename+'_ExtLCs')

	#### EXTRACT LIGHTCURVES (SELECTED)
	#ExtLCs(data['STAR'],data['RA'], data['DEC'], data['TIC'], 4, sspoc = data['SSPOC'], sffi = data['SFFI'],
	#		    thmask = data['MASK'], pdf_name=args.filename+'_ExtLCs')

	#### VISUALIZE LIGHTCURVES
	#VisLCs(star_ids=fdata['STAR'], plots_per_grid = 8, folder = 'EXTLC_R1', pdf_name = args.filename+'_LCs', 
	#								       pdf_save = True, inter = False)

	#### VISUALIZE PERIODOGRAMS
	#ProcessLevel(level='ls',star_ids=fdata['STAR'], rows_page = 6, cols_page = 2, folder = 'EXTLC_R1',
	#		        output_name = args.filename+'_LC', output_format = None, inter = False)

LS = Visualize(level='ls',star_ids=fdata['STAR'], rows_page = 5, cols_page = 5, 
			   output_name = args.filename+'_LS', coll_x = True, coll_y = True, output_format = None, inter = True)

ext_files = ['HAUCKE+19','ext_prop','FRASER+10']
	#ext_keys = ['TEFF','MDOT','VSINI','LOGQ','NABUN','LOGD','LOGLM']
	#ext_keys = ['LOGG','MASS','LOGL','VMAC','VMIC','VINF']
	#ext_keys = ['EW3995','EW4129','EW4131','EW4553','EW4568','EW4575']
	#for k in ext_keys:
	#	external = XmExtCol(inputf['RA'], inputf['DEC'],ext_files=ext_files,ext_col= k)
	#	fdata = hstack([fdata,external])

	#ext_keys = ['STAR'] + ['f_'+ k for k in ext_keys] + ext_keys


	#Correlation(ext_x = LS.rn_tab, ext_y = fdata[ext_keys], match_keys = 'STAR', 
	#		coll_x = True, coll_y = True, output_format = None, inter = True)



















