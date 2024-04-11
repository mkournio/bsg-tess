from constants import *
import os
import numpy as np
from functools import reduce
import itertools
import re
from astropy.table import Table

class FrequencyManager(object):

	def __init__(self, data, **kwargs):

		self.data = data
		self._freq_tab()	

		return 

	def _freq_tab(self):

	 	fr_tab = [['STAR','SECTORS','FF','A_FF','SNR_FF','FFR','A_FFR','HF']]

		# Collects, classifies all TESS stored frequencies

		for star in self.data['STAR']:

			star_files = [f for f in os.listdir(TESS_LS_PATH) if 'FREQ' in f and star in f]				
			
			st_sects = []
			st_ff = []
			st_ff_ampl = []
			st_ff_snr = []
			st_ffr = []
			st_ffr_ampl = []
			st_hf = []
			
			for sf in star_files:

				st_sects.append(sf.split('_')[1])

				sect_ff = []
				sect_ff_ampl = []
				sect_ff_snr = []
				sect_hf = []

				with open(TESS_LS_PATH+sf) as file_:
					first_line = file_.readline()
					for line in file_:
						fr, ampl, snr, ident = map(line.split().__getitem__,[1,2,5,6])
						if ident == '(F)':
							sect_ff.append(float(fr))
							sect_ff_ampl.append(float(ampl))
							sect_ff_snr.append(float(snr))
						elif ident == '(H)' or ident == '(HH)':
							sect_hf.append(float(fr))
				file_.close()

				st_ff.append(sect_ff)
				st_ff_ampl.append(sect_ff_ampl)
				st_ff_snr.append(sect_ff_snr)
				st_ffr.append([sect_ff[0]/sect_ff[1] if len(sect_ff) > 1 else np.nan])
				st_ffr_ampl.append([sect_ff_ampl[0]/sect_ff_ampl[1] if len(sect_ff) > 1 else np.nan])
				st_hf.append(sect_hf)

			fr_tab.append([star,st_sects,st_ff,st_ff_ampl,st_ff_snr,st_ffr,st_ffr_ampl,st_hf])

		self.freq_tab = Table(rows=fr_tab[1:], names=fr_tab[0])

		return 

class FrequencyIdentifier(object):

	# Given a list or array of frequencies, the function 'identify' returns the identification of frequencies
	# into fundamental, harmonics, and combinations of the fundamentals (currently combs of 2 first fundamental
	# are available). The frequency vector must be sorted by order of decreasing power amplitude, assuming that
	# the first frequency point in the vector is fundamental. Therefore, prior to calling the class
	# swapping between the first two points in order to start from the smallest value -typically taken as the fundamental- 
	# should be considered.

	def __init__(self, frequency_vector, resolution = 0.01):

		self.freq_vec = np.array(frequency_vector)
		self.resol = resolution	
		self.masks = [np.ones(len(frequency_vector),dtype=bool)]
		self.ident = np.copy(self.freq_vec)

	@property
	def ids(self):
		return self.identify()

	@property
	def ff(self):
		return self.freq_vec[self.ff_mask]

	@property
	def hf(self):
		return self.freq_vec[self.hf_mask]

	@property
	def cf(self):
		return self.freq_vec[self.cf_mask]

	@property
	def ff_mask(self):
		expr = r'^F\d{,3}$'
		return np.array([bool(re.search(expr,x)) for x in self.ids])	

	@property
	def hf_mask(self):
		expr = r'^\d{1,}F\d{,3}$'
		return np.array([bool(re.search(expr,x)) for x in self.ids])	

	@property
	def cf_mask(self):
		expr = r'^-?\d{1,}F\d{,3}\s-?\d{1,}F\d{,3}$'
		return np.array([bool(re.search(expr,x)) for x in self.ids])

	@property
	def indiv_mask(self):
		mask = np.zeros(len(self.ids), dtype = bool)
		u, indices = np.unique(self.ids, return_index=True)
		for i in indices: mask[i] = True
		return mask

	def identify(self):

		i = 0; fi = 0; fund = []
		while i < len(self.ident):

			comb = False				
			try:
				f = float(self.ident[i])
				self.cond = abs(self.freq_vec - f) < self.resol

				if len(fund) > 2 :
					bc = []			
					for c0, c1 in list(itertools.product(np.arange(-5,5,1), repeat=2)):
						if abs(f - c0 * fund[0] - c1 * fund[1]) < self.resol:
							bc.append((c0,c1))
							comb = True
					if comb:
						m_ind = np.argmin([abs(c[0]) + abs(c[1]) for c in bc])
						id_comb = '%sF0 %sF1' % (bc[m_ind][0],bc[m_ind][1])
						self._ident_update(id_comb)

				if not comb: 

					id_fund = 'F%d' % fi 
					self._ident_update(id_fund)

					self.cond = np.logical_or(self.freq_vec % f < self.resol, 
								  f - self.freq_vec % f < self.resol)
					mult =	np.array([str(int(round(x/f))) + id_fund for x in self.freq_vec])
					self._ident_update(mult)

					fund.append(f)
					fi += 1
			except:
				pass

			i += 1

		return self.ident

	def _ident_update(self, key_repl):	

			self.ident = np.where(self.cond & reduce(np.logical_and,self.masks), key_repl, self.ident)
			self.masks.append(~self.cond)
			
			return
