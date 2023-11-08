from plot_methods import GridTemplate
from functions import *
from astropy.table import join
import numpy as np

class corr_scatter(GridTemplate):

	def __init__(self, ext_x, ext_y, match_keys, **kwargs):

		self.ext_x = ext_x
		self.ext_y = ext_y
		self.match_keys = match_keys

		self.cols_x = [c for c in self.ext_x.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]
		self.cols_y = [c for c in self.ext_y.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]

		super(corr_scatter,self).__init__(rows_page = len(self.cols_x), cols_page = len(self.cols_y),
						 row_labels = self.cols_x, col_labels = self.cols_y, **kwargs)

		self._scatter_panel()
		self.GridClose()	

	def _scatter_panel(self):


		jtab = join(self.ext_x, self.ext_y, keys=self.match_keys)

		try:
			jtab['TEFF'] = np.log10(jtab['TEFF'])
		except:
			pass

		for c1 in self.cols_x :
			for c2 in self.cols_y :

				ax = self.GridAx()

				flag_fraser = jtab['f_'+c2] == 2
				ax.plot(jtab[c2][flag_fraser],jtab[c1][flag_fraser],'g^')

				flag_haucke = jtab['f_'+c2] == 1
				ax.plot(jtab[c2][flag_haucke],jtab[c1][flag_haucke],'ro')
				#ax.plot(jtab[c2],jtab[c1],'k.')				


class corr_hist(GridTemplate):

	def __init__(self, star_ids, ext_y, match_keys, type_data = 'frequency', **kwargs):

		self.star_ids = star_ids
		self.ext_y = ext_y
		self.match_keys = match_keys
		self.type_data = type_data	
		self.freqs = []
		self.snrs = []
		self.ident = []	

		self.cols_y = [c for c in self.ext_y.columns if not c.startswith(('STAR','e_','RA','DEC','f_'))]
		super(corr_hist,self).__init__(rows_page = 1, cols_page = 1,
						    **kwargs)

		self._get_data()

		self._hist_panel()
		self.GridClose()

	def _get_data(self):

		if self.type_data == 'frequency':	

			for star in self.star_ids:

				star_freqs = []
				star_snrs = []
				star_ident = []
				star_files = [f for f in os.listdir(TESS_LS_PATH) if 'FREQ' in f and star in f]				

				for sf in star_files:
					with open(TESS_LS_PATH+sf) as file_:
						file_.readline()
						for line in file_:
							fr, snr, ident = map(line.split().__getitem__,[1,5,6])
							star_freqs.append(fr)
							star_snrs.append(snr)
							star_ident.append(ident)
					file_.close()

				self.freqs.append(np.array(star_freqs).astype(np.float64))
				self.snrs.append(np.array(star_snrs).astype(np.float64))
				self.ident.append(np.array(star_ident))

		self.freqs = np.array(self.freqs)
		self.snrs = np.array(self.snrs)
		self.ident = np.array(self.ident)

		return

	def _hist_panel(self):
				ax = self.GridAx()		
				bins=np.arange(0.08,1.5,0.05);

				masked_freqs = mask_flatten(self.freqs,star_mask)

		
				ax.hist(masked_freqs, bins=bins, color= 'b')

	def _mask_flatten



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







