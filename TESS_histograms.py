import numpy as np
import os
import gc
from paths import *
from LPT_plotting import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

params = {'legend.fontsize': 10,
	 'font.size':  14,
         'axes.labelsize': 7,
         'axes.titlesize': 11,
         'xtick.labelsize': 12,
         'ytick.labelsize': 12}
plt.rcParams.update(params)

file_fold = os.listdir(TESS_FREQ_PATH)
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

