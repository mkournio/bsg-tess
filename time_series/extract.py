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
