from constants import *
import os
import numpy as np
from functools import reduce
import itertools
import re

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
