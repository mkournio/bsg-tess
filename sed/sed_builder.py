from plot_methods import GridTemplate
from astropy.table import Table, hstack
from scipy.interpolate import interp1d
from functions import *
from constants import *
import numpy as np
import os 

def synthetic_fluxes(synth_l, model_type = 'kurucz'):

	# Creates a tabular of synthetic SED fluxes over a vector of effective wavelengths. 
	# Model grids: Kurucz/PoWR

	model_tab = Table(names=['t','g']+synth_l)
	if model_type == 'kurucz':
		grid_teffs = sorted(set([float(f[1:6]) for f in os.listdir(KURUCZ_SED_PATH) if f.endswith('flx')]))
		grid_loggs = sorted(set([float(f[7:9]) for f in os.listdir(KURUCZ_SED_PATH) if f.endswith('flx')]))
	elif model_type == 'powr':
		grid_teffs = sorted(set([1e+3*float(f[5:7]) for f in os.listdir(POWR_SED_PATH) if f.startswith('ob')]))
		grid_loggs = sorted(set([float(f[8:10]) for f in os.listdir(POWR_SED_PATH) if f.startswith('ob')]))

	for t in grid_teffs:
		for g in grid_loggs:
			try:
					sed = sed_read(t,g,model_type,synth_l)
					model_tab.add_row(np.append([t,g],sed))				
			except:
				pass

	return model_tab


def sed_read(temp,logg,models,wave):

	# Read SED files and interpolate flux at given list of wavelengths. 

	if models == 'kurucz':
		sed_name = 't%05dg%02d.flx' % (temp,logg)
  		tab = np.genfromtxt(KURUCZ_SED_PATH+sed_name, names = "wave, flux")
	elif models == 'powr':
		sed_name = 'ob-i_%2d-%2d_sed.txt' % (1e-3*temp,logg)
  		tab = np.genfromtxt(POWR_SED_PATH+sed_name, names = "wave, flux")
		tab['wave'] = 10 ** tab['wave']
		# Conversion into surface flux - raw is flux received at D = 10 pc
		for l in open(POWR_SED_PATH+'modelparameters.txt'):
			if '%2d-%2d'% (1e-3*temp,logg) in l :
				lstar = 10**float(l.split()[5])
				rad2 = lstar * ((float(l.split()[1])/5778.)**-4)
				scalar = rad2 * (1e-1 *  RSUN_TO_PC)**2
		tab['flux'] = (10 **  tab['flux']) * (scalar**-1)
		
	f = interp1d(tab['wave'],tab['flux'],kind='linear')

 	return f(wave)

'''
class SEDBuilder(GridTemplate):

	def __init__(self, photo_tab, filt_dict, fit_dict = None, **kwargs):

		super(SEDBuilder,self).__init__(fig_xlabel=PLOT_XLABEL['sed'],
					        fig_ylabel=PLOT_YLABEL['sed'],**kwargs)
		self.photo_tab = photo_tab
		self.filt_dict = filt_dict

		if fit_dict != None :

			fit_keys = ['MODELS','FIT_BODIES','ST_PROP']

			if not all([x in fit_dict for x in fit_keys]):
				raise ValueError("One or more of the fit columns are missing.")
			self.__dict__.update(fit_dict)

			VALIDATE HERE

			from lmfit import minimize, Parameters, Parameter
			from bisect import bisect

			self.minimize, self.Parameter, self.Parameters = minimize, Parameter, Parameters
			self.bisect = bisect

		self._sed_plot()

		self.GridClose()

		return

	def _sed_plot_single(self, t_ind, star_label):

		ax = self.GridAx()

		if self.fit_dict != None : self.fit_l=[]; self.fit_f=[]; self.fit_ef=[]

 		for k, v in self.filt_dict.items() :

			photo_star = []; photo_star_err = []
			for c, ec in zip(v['f'],v['e']):

				if c in self.photo_tab.columns:
					photo_star.append(self.photo_tab[c][t_ind])
				elif c+'/'+k in self.photo_tab.columns:
					photo_star.append(self.photo_tab[c+'/'+k][t_ind])	
				else:
					photo_star.append(np.ma.masked)

				if ec in self.photo_tab.columns:
					photo_star_err.append(self.photo_tab[ec][t_ind])
				elif ec+'/'+k in self.photo_tab.columns:
					photo_star_err.append(self.photo_tab[ec+'/'+k][t_ind])	
				else:
					photo_star_err.append(np.ma.masked)				

			zeros = v['z']
			lambdas = v['l']

			photo_star_flx, photo_star_eflx = mg2flux(photo_star,photo_star_err,zeros,lambdas,k)

			if not np.ma.array(photo_star).all() is np.ma.masked:
				print k, np.array([lambdas, photo_star, photo_star_err, photo_star_flx, photo_star_eflx]).T
				ax.plot(lambdas,photo_star_flx,v['mrk']['m'],c=v['mrk']['c'],markersize=v['mrk']['s'],label=k)
  				ax.errorbar(lambdas,photo_star_flx,photo_star_eflx,ecolor=v['mrk']['c'],elinewidth=1, capsize=0, ls='none')

			if self.fit_dict != None and v['fit'] == 1:
				self.fit_l.extend(np.array(lambdas))
				self.fit_f.extend(np.array(photo_star_flx))
				self.fit_ef.extend(np.array(photo_star_eflx))
				#ax.scatter(np.array(lambdas)[fit_ind],3*np.array(photo_star_flx)[fit_ind],marker=u'$\u2193$',\
				#	color='c',s=17)

		if self.fit_dict != None: 
			try:
				self.teff = self.teff_tab[t_ind][0]
				self.dist = self.dist_tab[t_ind][0]
				self.logg = int(10*self.logg_tab[t_ind][0])

				if self.teff > 15000. :
					self.models_key = 'powr'
				else:
					self.models_key = 'kurucz'

				self.mod_tab = self.mod_tabs[self.models_key]
				self.mod_tarr = sorted(set(self.mod_tab['t']))

				fit_prop = self._fit()
				self._plot_model(ax,fit_prop)

			except:
				pass

		ax.legend(loc='upper right')
		ax.set_xscale('log'); ax.set_yscale('log')
		ax.set_title(star_label)

		return

	def _fit(self):

		params = self.Parameters()
		for fb in self.fit_bodies :

			self.fit_body = fb

			bod_prop = SED_BODIES[fb]
			for p in bod_prop:
				params[p] = self.Parameter(**bod_prop[p])

			if fb == 'psph' :
				params['teff'].set(value = self.teff)
				params['dist'].set(value = self.dist)
				if len(self.fit_bodies) == 1 :
					self.mask = map(lambda x: x < LAMBDA_FIT_THRES, self.fit_l)
				else:
					self.mask = map(lambda x: x < min(1.5e+4,LAMBDA_FIT_THRES), self.fit_l)

			elif fb == 'hd' and 'wd' in self.fit_bodies:
				self.mask = map(lambda x: x < min(3e+4,LAMBDA_FIT_THRES), self.fit_l)

			else:
				self.mask = map(lambda x: x <= LAMBDA_FIT_THRES, self.fit_l)

			result = self.minimize(self._resid, params, method = 'least_squares')

			for k in params:
				val=result.params[k].value
				params[k].set(value = val)
				params[k].set(min = 0.8 * val)
				params[k].set(max = 1.2 * val)
 		
		return result

	def _resid(self, params):

		tr_l = np.array(self.fit_l,dtype=object)[self.mask]
		tr_f = np.array(self.fit_f)[self.mask]
		tr_l_mu = tr_l / 10000.

  		rad = params['rad'].value
		ext = params['ext'].value
		dist = params['dist'].value
		teff = params['teff'].value
		
		t_low, t_up, q, logg_gr = self._get_tvals(teff,self.logg)	
		model_low = self.mod_tab[(self.mod_tab['t'] == t_low) & (self.mod_tab['g'] == logg_gr)]
		model_up = self.mod_tab[(self.mod_tab['t'] == t_up) & (self.mod_tab['g'] == logg_gr)]

		synth_f = []
		for l in tr_l:
			synth_flux = (1-q)*model_low[str(l)] + q*model_up[str(l)]
			synth_f.append(synth_flux[0])
		synth_f = np.array(synth_f)

		scalar = self._get_scalar(rad, dist)
		scaled_f = scalar * synth_f * (10**(-0.4*red_coeff(tr_l_mu,ext,3.1))) 

		try:
			hds = params['hds'].value
			hdt = params['hdt'].value
			hd_f = (IR_SCALAR * hds) * planck(tr_l,hdt,kind='grey')
			scaled_f += hd_f
	
			wds = params['wds'].value
			wdt = params['wdt'].value
			wd_f = (IR_SCALAR * wds) * planck(tr_l,wdt,kind='grey')
			scaled_f = scaled_f + wd_f	
		except:
			pass

		resid = np.true_divide(abs(scaled_f-tr_f),tr_f)
		resid = resid[~np.isnan(resid)]

		return resid

	def _get_tvals(self,teff,logg):

		t_ind = self.bisect(self.mod_tarr, teff)
		t_low = self.mod_tarr[t_ind-1]
		t_up = self.mod_tarr[t_ind]
		q =  (teff - t_low) / (t_up - t_low)

		logg_gr = logg
		cond = False
		while not cond:
			cond1 = (self.mod_tab['t'] == t_low) & (self.mod_tab['g'] == logg_gr)
			cond2 = (self.mod_tab['t'] == t_up) & (self.mod_tab['g'] == logg_gr)
			cond = (True in cond1) and (True in cond2)
			logg_gr += 1

		return int(t_low), int(t_up), q, logg_gr-1		

	def _get_scalar(self,radius, dist):
		return (radius * RSUNM / (float(dist) * PC_TO_M))**2

	def _plot_model(self,axis,fit_prop):

		rad = float(fit_prop.params['rad'].value) 
		ext = float(fit_prop.params['ext'].value)
		dist = float(fit_prop.params['dist'].value)
		teff = float(fit_prop.params['teff'].value)
		logL = np.log10( (rad**2) * ((teff/5778.)**4) )

		#print fit_prop.params
		
		w = np.append(np.arange(1000,100000,1),[200000,400000])
		t_low, t_up, q, logg_gr = self._get_tvals(teff,self.logg)
		sed_low = sed_read(t_low,logg_gr,self.models_key,w)
		sed_up = sed_read(t_up,logg_gr,self.models_key,w)

		sed_interp = (1-q) * sed_low + q * sed_up

		scalar = self._get_scalar(rad, dist)
		w_mu = w/10000.
		scaled_sed = scalar * sed_interp * (10**(-0.4*red_coeff(w_mu,ext,3.1)))
		scaled_sed_unr = scalar * sed_interp 
		axis.plot(w,scaled_sed_unr,'k:')

		hdtext, wdtext = '',''
		try:
			hdt = float(fit_prop.params['hdt'].value)
			hds = float(fit_prop.params['hds'].value)
			hd_f = (IR_SCALAR * hds) * planck(w,hdt,kind='grey')
			scaled_sed += hd_f
			hd_f[w < 5e+3] = np.nan			
			axis.plot(w,hd_f,'r')
			hdtext = 'T$_\mathrm{hd}$ [K] = %s\n' % roundup(hdt)

			wdt = float(fit_prop.params['wdt'].value)
			wds = float(fit_prop.params['wds'].value)
			wd_f = (IR_SCALAR * wds) * planck(w,wdt,kind='grey')
			scaled_sed += wd_f
			wd_f[w < 1e+4] = np.nan			
			axis.plot(w,wd_f,'b')	
			wdtext = 'T$_\mathrm{wd}$ [K] = %s\n' % roundup(wdt)
		except:
			pass

		axis.plot(w,scaled_sed,'k')
	

		vtext = '%s\n' % MODEL_LEG[self.models_key] + \
			'T$_\mathrm{eff}$ [K] = ' + self._pr(int(teff),True) + \
			'logg [dex] = ' + self._pr(1e-1*logg_gr,False,note='(sp %s)' % (1e-1*self.logg)) + \
			'R [$R_{\odot}$] = ' + self._pr(int(round(rad)),False) + \
			'D [pc] = ' + self._pr(int(dist),True) + \
			'log(L/$L_{\odot}$) = ' + self._pr(round(logL,2),False) + \
			'A$_{V}$ [mag] = ' + self._pr(round(ext,2),False)
	
		axis.text(0.02,0.005,vtext,size=10,transform=axis.transAxes)
		axis.text(0.7,0.1,hdtext+wdtext,size=10,transform=axis.transAxes)
		
		return

	def _pr(self,value,fixed,note=''):
		if fixed:
			return '%s [F]\n' % value 
		else:	
			return '%s %s\n' % (value,note)
			
'''
