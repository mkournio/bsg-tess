from constants import *
import os
import numpy as np
import numpy.ma as ma
from functools import reduce

class FreqIdent(object):

	# Class for the identification of frequencies of a star object

	def __init__(self,star_id):

		self.star_id = star_id

	def _read_freq(self):

		star_files = [f for f in os.listdir(TESS_LS_PATH) if 'FREQ' in f and self.star_id in f]				

		star_freqs = []
		star_snrs = []
		star_ident = []

		im = np.arange(2,16,1)

		for sf in star_files:

			sect_f = []
			sect_snr = []

			with open(TESS_LS_PATH+sf) as file_:
				first_line = file_.readline()
				rayleigh = float(first_line.split()[1])
				for line in file_:
					fr, snr, ident = map(line.split().__getitem__,[1,5,6])
					star_freqs.append(float(fr))

					sect_f.append(float(fr))
					sect_snr.append(float(snr))
					
					star_snrs.append(float(snr))
					star_ident.append(ident)
			file_.close()

			sect_f = np.array(sect_f)
			sect_snr = np.array(sect_snr)

			sect_f = sect_f[ sect_f > 2 * rayleigh]

			m_ind = np.argmin(sect_f[:2])
			sect_f[0], sect_f[m_ind] = sect_f[m_ind], sect_f[0]
			sect_snr[0], sect_snr[m_ind] = sect_snr[m_ind], sect_snr[0]

			ident = np.copy(sect_f)

			i = 0; fi = 0

			print ident
			masks = [np.ones(len(sect_f),dtype=bool)]

			while i < len(ident):
				
				try:
					f = float(ident[i])

					idtx = 'F%d' % fi

					cond = abs(sect_f - f) < rayleigh 
					ident = np.where(cond & reduce(np.logical_and,masks), idtx, ident)
					#print cond & reduce(np.logical_and,masks), ident

					masks.append(~cond)

					mult =	np.array([str(int(round(x/f))) + idtx for x in sect_f])
					cond = np.logical_or(sect_f % f < rayleigh, f - sect_f % f < rayleigh)
					ident = np.where(cond & reduce(np.logical_and,masks), mult, ident)
					#print cond & reduce(np.logical_and,masks), ident

					masks.append(~cond)


					fi += 1
				except 
					pass

				i += 1

			print ident


	def _flconv(self,x):
		try:
			return float(x)
		except:
			return x
			

		#return star_freqs, star_snrs, star_ident



def find_harmon(v,e):

	'''
	Calculates the harmonics/combinations from a vector of frequencies
	'''
	v = np.array(v)
	ident_v = np.zeros(len(v),dtype=np.dtype('U100'))

	f = np.where(v > 2 * e)[0]
	t = np.where(v <= 2 * e)[0]

	np.put(ident_v,t,'(BTH)')

	ratio = v[f[0]]/v[f[1]]
	error = ratio * e * np.sqrt( (1/v[f[0]]**2) + (1/v[f[1]]**2) )	

	if (abs(ratio - round(ratio)) < error) and (round(ratio) == 1) :
			id0 = '(F)'
			id1 = '(U) F%d' % f[0]
	elif abs(ratio - round(ratio)) < error :
			id0 = '(H) %dF%d' % (round(ratio),f[1])
			id1 = '(F)'
	elif abs((1/ratio) - round(1/ratio)) < error * (1/ratio**2) :
			id0 = '(F)'
			id1 = '(H) %dF%d' % (round(1/ratio),f[0])
	else:
			id0 = '(F)'
			id1 = '(F)'

	ident_v[f[0]] = id0
	ident_v[f[1]] = id1

	for k in f[2:] :

		ind_k = np.where(f==k)[0][0]

		minim = 100.
 		oc1 = '0F'; oc2 = '0F'
 		for i in f[:ind_k-1] :
			for j in f[i+1:ind_k] :
					for c1 in np.arange(-5,10,1):
    						for c2 in np.arange(-5,10,1):

				        		if (abs(c1*v[i] + c2*v[j] - v[k]) < np.sqrt(2) * e) and \
								(abs(c1) + abs(c2) < minim): 

								minim = abs(c1) + abs(c2)
								oc1 = '%sF%s' % (c1,i); oc2 = '%sF%s' % (c2,j)

	 	sum_coeff = abs(float(oc1.split('F')[0]))+abs(float(oc2.split('F')[0]))
		oc1 = trans_freq_text(oc1)
		oc2 = trans_freq_text(oc2)

	 	text = '%s %s' % (oc1,oc2)

	 	if sum_coeff == 0:
			text = '(F)'
	 	elif sum_coeff == 1:
			text = '(U) ' + text
	 	elif sum_coeff > 1 and ('' in [oc1,oc2]):
			text = '(H) ' + text
	 	else:
			text = '(C) ' + text

	 	ident_v[k] = text 

	return ident_v




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
