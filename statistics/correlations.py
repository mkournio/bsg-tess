from plot_methods import GridTemplate
from functions import *
from constants import *
from astropy.table import join
import numpy as np
from scipy.stats import norm, gamma 
from scipy.stats.mstats import pearsonr, spearmanr, describe

class corr_scatter(GridTemplate):

	def __init__(self, ext_x, ext_y, match_keys = 'STAR', auto = False, **kwargs):

		self.ext_x = ext_x
		self.ext_y = ext_y
		self.match_keys = match_keys
		self.auto = auto

		self.cols_x = [c for c in self.ext_x.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]
		self.cols_y = [c for c in self.ext_y.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]

		super(corr_scatter,self).__init__(rows_page = len(self.cols_x), cols_page = len(self.cols_y),
						 row_labels = self.cols_x, col_labels = self.cols_y, **kwargs)

		self._scatter_panel()
		self.GridClose()	

	def _scatter_panel(self):

		if self.auto:
			jtab = self.ext_x
		else:
			jtab = join(self.ext_x, self.ext_y, keys=self.match_keys)

		if 'TEFF' in jtab.columns: jtab['TEFF'] = np.log10(jtab['TEFF'])

		for c1 in self.cols_x :

			nx = mask_outliers(jtab[c1], m = 4)

			for c2 in self.cols_y :

				ax = self.GridAx()

				flag_fraser = jtab['f_'+c2] == 2
				ax.plot(jtab[c2][flag_fraser],jtab[c1][flag_fraser],'^',color='limegreen')

				flag_haucke = jtab['f_'+c2] == 1
				ax.plot(jtab[c2][flag_haucke],jtab[c1][flag_haucke],'ro')				

				ny = mask_outliers(jtab[c2], m = 4)

				ax.plot(ny,nx,'k.',markersize=2)

				if not (nx.mask.all() or ny.mask.all()):
					pears_val = pearsonr(x=nx,y=ny)[0]
					spear_val = spearmanr(x=nx,y=ny)[0]
					ax.text(0.05,0.30,'%.2f' % pears_val,size = 6,color='c',transform=ax.transAxes)
					ax.text(0.05,0.15,'%.2f' % spear_val,size = 6,color='b',transform=ax.transAxes)

				#ax.plot(jtab[c2],jtab[c1],'k.')	
			
class autocorr_scatter(corr_scatter):

	def __init__(self, vector, **kwargs):

		super(autocorr_scatter,self).__init__(ext_x = vector, ext_y = vector, auto = True, **kwargs)
		return


class corr_hist(GridTemplate):

	def __init__(self, ext_tab, snr_show = 5, type_data = 'frequency', **kwargs):

		self.ext_tab = ext_tab
		self.type_data = type_data
		self.snr_show = snr_show	
		self.freqs = []
		self.snrs = []
		self.ident = []	

		self.cols_y = [c for c in self.ext_tab.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]
		super(corr_hist,self).__init__(rows_page = len(self.cols_y), cols_page = 3, row_labels = self.cols_y, 
					       fig_xlabel = PLOT_XLABEL['ls'], 
						sup_xlabels = ['Fundamental (S/N > %s)' % self.snr_show,'Fundamental', 'Harmonics'],
					       **kwargs)

		self._get_data()

		self._hist_panel()
		self.GridClose()

	def _get_data(self):

		if self.type_data == 'frequency':	

			for star in self.ext_tab['STAR']:

				star_freqs = []
				star_snrs = []
				star_ident = []
				star_files = [f for f in os.listdir(TESS_LS_PATH) if 'FREQ' in f and star in f]				

				for sf in star_files:
					with open(TESS_LS_PATH+sf) as file_:
						file_.readline()
						for line in file_:
							fr, snr, ident = map(line.split().__getitem__,[1,5,6])
							star_freqs.append(float(fr))
							star_snrs.append(float(snr))
							star_ident.append(ident)
					file_.close()

				self.freqs.append(star_freqs)
				self.snrs.append(star_snrs)
				self.ident.append(star_ident)

			self.freqs = fill_array(self.freqs)
			self.snrs = fill_array(self.snrs)
			self.ident = fill_array(self.ident)

			return

	def _hist_panel(self):
						
		bins=np.arange(0.05,1.15,0.04);


		for prop in self.cols_y:

			ax_hsn = self.GridAx()
			ax_all = self.GridAx()
			ax_depf = self.GridAx()

			nanprop = self.ext_tab[prop].filled(np.nan)

			brp_ini = [np.nanpercentile(nanprop, 33.33),np.nanpercentile(nanprop, 66.66)] 

			tb = [-np.inf] + sorted(brp_ini) + [+np.inf]

			for i in range(len(tb)-1):

				mask_star = np.logical_and(tb[i] < self.ext_tab[prop],self.ext_tab[prop] < tb[i+1])
				mask_star = np.logical_and(mask_star,self.ext_tab[prop].mask == False)

				bin_freqs = self.freqs[mask_star]
				bin_snrs = self.snrs[mask_star]
				bin_ident = self.ident[mask_star]

				for s, ax in zip([self.snr_show,0],[ax_hsn,ax_all]):
							
					mask_fund = np.logical_and(bin_snrs > s, bin_ident == '(F)')
					hdata = bin_freqs[mask_fund]

					hdata_fit = [x for x in hdata if abs(x - np.mean(hdata)) < 3 * np.std(hdata)]
					#(mu, sigma) = norm.fit(hdata_fit)
					gparam = gamma.fit(hdata_fit, floc=0)

					if len(hdata) > 0: 
						values, _, _ =ax.hist(hdata, histtype='step', lw = 2, color = HIST_COLORS[i], 
							bins=bins, label = self._hist_legend(i,STY_LB[prop],tb[i],tb[i+1],FMT_PR[prop]))
								
						xfit = np.linspace(min(hdata), max(hdata), 100)
						area = sum(np.diff(bins) * values)

						if s > 0 : print prop, i, area

						#ax.plot(xfit,norm.pdf(xfit,mu,sigma)*area, lw=2.5, color = HIST_COLORS[i])
						ax.plot(xfit,gamma.pdf(xfit,*gparam)*area, lw=2.5, color = HIST_COLORS[i])

						ax.set_ybound(lower=0.1)
						ax.set_xbound(lower=0.01)		


				mask_depf = np.logical_and(bin_snrs > 0, bin_ident == '(H)')
				data_depf = bin_freqs[mask_depf]
				ax_depf.hist(data_depf, alpha = 0.4, histtype='step', lw = 2, color = HIST_COLORS[i],bins=bins)

			ax_hsn.text(0.7,0.3,'stbin_size %s' % mask_star.sum(), size = 8, transform=ax_hsn.transAxes)
			ax_all.legend(loc="upper right")	


	def _hist_legend(self,i,prop,low,high,frm):

		if i == 0:
			return '%s < ' % prop + frm % high
		elif i == 1:
			return frm % low + ' < %s < ' % prop + frm % high
		elif i == 2:
			return '%s > ' % prop + frm % low


'''

			

			


	def _corr_panel(self):



data = ascii.read('LaPlata_input')
data['RA'], data['DEC'] = to_deg(data['RA'],data['DEC'])

def freq_type(data,t):
	
	fund=[]; harmon=[]; comb=[]; unres=[]

	for star, ident in zip(data['STAR'],data['TCATG']):

  	   if ident == t:	
		star_files = [j for j in file_fold if (star in j) and ('FREQ' in j) and ('PL' not in j)]
		if len(star_files) > 0 :
			
			s_fund=[]; s_harmon=[]; s_comb=[]; s_unres=[]; 

			for s in star_files :
				with open(TESS_FREQ_PATH+s) as f:

					for line in f:

						freq = float(line.split()[1])

						if '(F)' in line:
							s_fund.append(freq)
							#if freq > 0.43 and ident == 'eB' : print star
						elif '(H)' in line:
							s_harmon.append(freq)
						elif '(C)' in line:
							s_comb.append(freq)
						elif '(U)' in line:
							s_unres.append(freq)
				f.close()

			for k in set(np.digitize(s_fund,bins)):
				mask = np.digitize(s_fund,bins) == k
				fund.append(np.mean([x for x,y in zip(s_fund,mask) if y]))
			for k in set(np.digitize(s_harmon,bins)):
				mask = np.digitize(s_harmon,bins) == k
				harmon.append(np.mean([x for x,y in zip(s_harmon,mask) if y]))
			for k in set(np.digitize(s_comb,bins)):
				mask = np.digitize(s_comb,bins) == k
				comb.append(np.mean([x for x,y in zip(s_comb,mask) if y]))

	return fund, harmon, comb

def fit_function(x, B, mu, sigma):
    return (B * np.exp(-1.0 * (x - mu)**2 / (2 * sigma**2)))

bins=np.arange(0.08,5,0.05);
binscenters = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])

efund, eharmon, ecomb = freq_type(data,'eB')
mfund, mharmon, mcomb = freq_type(data,'mB')
lfund, lharmon, lcomb = freq_type(data,'lB')

fig,[ax1,ax2,ax3] = plt.subplots(3,1,sharex=True)

ax1.set_title('late O/early B')
ax2.set_title('mid-B')
ax3.set_title('late B/early A')

ax1.hist(efund, bins=bins, alpha = 1.0, color= 'b',label='Fund. freq.')
ax2.hist(mfund, bins=bins, alpha = 1.0, color= 'g')
ax3.hist(lfund, bins=bins, alpha = 1.0, color= 'r')

edata, _, _ = ax1.hist(eharmon, bins=bins, histtype=u'step', color= 'b',label='Harmonics')
mdata, _, _ = ax2.hist(mharmon, bins=bins, histtype=u'step', color= 'g')
ldata, _, _ = ax3.hist(lharmon, bins=bins, histtype=u'step', color= 'r')
e_popt, _ = curve_fit(fit_function, xdata=binscenters, ydata=edata, p0=[2.0, 0.5, 1.0]); print e_popt
m_popt, _ = curve_fit(fit_function, xdata=binscenters, ydata=mdata, p0=[2.0, 0.5, 1.0]); print m_popt
l_popt, _ = curve_fit(fit_function, xdata=binscenters, ydata=ldata, p0=[2.0, 0.5, 1.0]); print l_popt
ax1.plot(bins, fit_function(bins, *e_popt), 'b')
ax2.plot(bins, fit_function(bins, *m_popt), 'g') 
#ax2.plot(bins/8, fit_function(bins, *m_popt), 'w--'); 
ax3.plot(bins, fit_function(bins, *l_popt), 'r')
#ax3.plot(bins/7, fit_function(bins, *l_popt), 'w--'); 

ax1.hist(ecomb, bins=bins, alpha = 0.2, color= 'b',label='Comb. freq.')
ax2.hist(mcomb, bins=bins, alpha = 0.2, color= 'g')
ax3.hist(lcomb, bins=bins, alpha = 0.2, color= 'r')

plt.xlim(bins[0]-0.01,bins[-1]+0.01)
ax1.legend(loc=1)

fig.text(0.5,0.06,'Frequency (c/d)', ha="center", va="center")
fig.text(0.06,0.5,'Counts', ha="center", va="center", rotation=90)

plt.show()


'''







