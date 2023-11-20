from plot_methods import GridTemplate
from functions import *
import lightkurve as lk
import numpy as np
import os 


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

