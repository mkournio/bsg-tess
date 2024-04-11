from plot_methods import GridTemplate
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from functions import *
import numpy as np
import lightkurve as lk
import os
import pickle
from scipy import stats

class Processing(GridTemplate):
	
	def __init__(self, data, level, load_rn_pickle = False, load_mt_pickle = False, **kwargs):

		if load_rn_pickle and level == 'ls' :

			self.rn_tab = pickle.load(open(PICKLE_PATH+'rn.pkl','rb'))
			print 'Loaded pickle: red noise properties - no plot is generated'

			return

		elif load_mt_pickle and level == 'lc' :

			self.mt_tab = pickle.load(open(PICKLE_PATH+'mt.pkl','rb'))
			print 'Loaded pickle: light curve metrics - no plot is generated'

			return

		else:
			self._validate()
			self.data = data
			self.level = level

			if level == 'lc':
				self.mt_tab = Table(names=('STAR','MAD','e_MAD', 'SVAR','e_SVAR','ZCROSS','e_ZCROSS','PSI','e_PSI','ETA','e_ETA','SKEW','e_SKEW'), dtype=('<U16','f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4', 'f4'))

			if level == 'ls': 
				kwargs['coll_x'] = True
				kwargs['coll_y'] = True
		 		self.rn_tab = Table(names=('STAR','w','zero','tau','gamma','e_w','e_zero','e_tau','e_gamma'), 
				    dtype=('<U16','f4','f4', 'f4', 'f4','f4', 'f4','f4', 'f4'))

			super(Processing,self).__init__(params = PLOT_PARAMS[level], fig_xlabel=PLOT_XLABEL[self.level],
					       fig_ylabel=PLOT_YLABEL[self.level],**kwargs)
			self._process()
			self.GridClose()
			if level == 'ls' and kwargs['output_format'] == None: 
				pickle.dump(self.rn_tab,open(PICKLE_PATH+'rn.pkl','wb'))
				print 'Saved pickle: red noise properties'	

			if level == 'lc' and kwargs['output_format'] == None: 
				pickle.dump(self.mt_tab,open(PICKLE_PATH+'mt.pkl','wb'))
				print 'Saved pickle: light curve metrics'

			return

	def _validate(self):
		pass

	def _process(self):

	        for star in self.data['STAR'] :
			if self.level == 'lc':
				lc_files = [j for j in os.listdir(TESS_LC_PATH) if (star in j)]
				mt = self._process_lc(star,lc_files)
				self.mt_tab.add_row(np.append(star,mt))

			elif self.level == 'ls':				
				rnopt = self._process_ls(star)
				self.rn_tab.add_row(np.append(star,rnopt))			
		return	

	def _process_ls(self, star):

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
			ax.plot(ls_v, self._redn(ls_v,*rnopt),'k-')

			star_mask = self.data['STAR'] == star
			vcrit = self.data['VCRIT'][star_mask]
			vsini = self.data['VSINI'][star_mask]
			radius = self.data['S_RAD'][star_mask]

			fund, harm, comb = self._get_freqs(TESS_LS_PATH + freq_files[np.argmin(rn_chi2s)])
			for l in fund : ax.axvline(l,color='b',ls='-',lw=2.5,ymin=0.80,ymax=1)
			for l in harm : ax.axvline(l,color='g',ls='--',lw=2,ymin=0.84,ymax=1)
			for l in comb : ax.axvline(l,color='k',ls=':',lw=1.5,ymin=0.87,ymax=1)

			rconst = 1. / (2 * PI * RSUNKM * S_TO_D)
			ax.axvspan(rconst * vsini / radius, rconst * vcrit / radius, color='0.85')
			try: 
			  tx = '$v$sin$i$ = %d kms$^{-1}$\n$R$ = %d R$_{\odot}$' % (vsini,radius)
			except:
			  tx = ''
			ax.text(0.65,0.65,tx,color='k',size=9, transform=ax.transAxes)

			ax.set_ylim(7e-7,2e-1) ;  ax.set_yscale('log')  
			ax.set_xlim(2. / TESS_WINDOW_D,25); ax.set_xscale('log') 

			ax.text(0.03,0.05,'%s (%s)' % (star,vis_sect),color='k',transform=ax.transAxes)

			rnerr[0] = rnerr[0]/(np.log(10)*rnopt[0])
			rnerr[1] = rnerr[1]/(np.log(10)*rnopt[1])
			rnopt[:2] = np.ma.log10(rnopt[:2])

			return rnopt, rnerr

	def _process_lc(self, star, lc_files):

			ax = self.GridAx()
			divider = make_axes_locatable(ax)

			ymin, ymax = self._lcs_glob_minmax(lc_files)
			ym = 1.1 * max(abs(ymax),abs(ymin))

			g_lc_files = group_consec_lcs(lc_files)	
			len_gls = len(g_lc_files)

			k = 0
			mad = []; zcross = []; psi = []; eta = []; svar = []; skew = []
			for gl in g_lc_files:

				g_times = np.empty(0); g_fluxes = np.empty(0); g_sects = []
				for l in gl:
					time, flux = lc_read(TESS_LC_PATH + l)
					offset = np.nanmedian(flux)
					g_sects.append(sect_from_lc(l))	
					g_times = np.append(g_times,time)
					g_fluxes = np.append(g_fluxes,flux - offset)				

				g_sect_tx = '+'.join([str(x) for x in g_sects])

				ax.plot(g_times, g_fluxes,'k'); ax.invert_yaxis(); ax.set_ylim(ym,-ym)
				ax.xaxis.set_tick_params(direction="in")
				#ax.yaxis.set_tick_params(labelsize=SIZE_XLABEL_SUB)
				if k == 0: ax.text(0.05,0.85,star,color='r',fontsize=SIZE_FONT_SUB,transform=ax.transAxes)
				ax.text(0.4,0.05,g_sect_tx,color='b',fontsize=SIZE_FONT_SUB,transform=ax.transAxes)

				x1_p = (2*g_times[0]+g_times[-1])/3.
				x2_p = (2*g_times[-1]+g_times[0])/3.
				ax.set_xticks([round_to(x1_p,5)[0], round_to(x2_p,5)[1]])

				if (len_gls > 1) and (k < len_gls - 1):
					d=.02
					kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
					ax.plot((1,1),(-d,+d), **kwargs) 
					ax.plot((1,1),(1-d,1+d), **kwargs) 
					ax.spines['right'].set_visible(False)

					ax = divider.append_axes("right", size="100%", pad=0.14)

					kwargs.update(transform=ax.transAxes) 
					ax.plot((0,0),(-d,+d), **kwargs) 
					ax.plot((0,0),(1-d,1+d), **kwargs) 
					ax.spines['left'].set_visible(False)
					ax.tick_params(labelleft = False) 

				nan_mask = np.isnan(g_fluxes)
				g_fluxes = g_fluxes[~nan_mask]
				
				mad.append( np.median(np.absolute(g_fluxes - np.median(g_fluxes))) )
				svar.append( np.std(g_fluxes) )

				zcross.append( k_crossing(g_times, g_fluxes, 1)[0] )
				psi.append( psi_sq(g_times, g_fluxes, 5) )
                                eta.append( np.ma.log10(eta_e(g_times, g_fluxes)) )	
				skew.append( stats.skew(g_fluxes) )

				k += 1

			return [np.mean(mad), np.std(mad), np.mean(svar), np.std(svar), np.mean(zcross), np.std(zcross),  np.mean(psi), np.std(psi), np.mean(eta), np.std(eta), np.mean(skew), np.std(skew)]

	def _get_chi2(self,path_freq):

		with open(path_freq) as f:
    			for line in f:
        			pass
    		
		return line.split()[3]

	def _get_freqs(self,path_freq):

		fund = []
		harm = []
 		comb = []
		with open(path_freq) as f:
    			first_line = f.readline()
    			for line in f:
				freq = float(line.split()[1])
				if '(F)' in line :
					fund.append(freq)
				elif ('(H)' in line) or ('(HH)' in line):
					harm.append(freq)
				elif '(C)' in line:
					comb.append(freq)

		return fund, harm, comb		
		

	def _redn(self,x,w,zero,tau,gamma):
		x = np.array(x)
		return w + ( zero / ( 1 + (2 * 3.14 * tau * x)**gamma))

	def _lcs_glob_minmax(self,lc_files):

		minfs=[]; maxfs=[]	
		for l in lc_files:
			time, flux = np.loadtxt(TESS_LC_PATH + l, delimiter=' ', usecols=(0, 1), unpack=True)
			time_dif = time[1:] - time[:-1]; ind_nan = np.argmax(time_dif)
			flux[ind_nan-3:ind_nan+15] = np.nan; flux[:5] = np.nan; flux[-30:] = np.nan
			minfs.append(np.nanmin(flux)); maxfs.append(np.nanmax(flux))

		return min(minfs), max(maxfs)





