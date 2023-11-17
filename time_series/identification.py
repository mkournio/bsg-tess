from constants import *
import os
import numpy as np
import numpy.ma as ma
from functools import reduce
import itertools

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

			print sect_f, FrequencyIdentify(sect_f, resolution = rayleigh).identify()



			return 1,2,3


class FrequencyIdentify(object):

	def __init__(self, frequency_vector, resolution = 0.01):

		self.freq_vec = np.array(frequency_vector)
		self.resol = resolution	
		self.ident = np.copy(frequency_vector)
		self.masks = [np.ones(len(frequency_vector),dtype=bool)]		

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

