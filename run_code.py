#%%
import args
from astropy.table import Table
from astropy.io import ascii
from functions import *
from time_series import *
from statistics import *
from sed import *
from constants import *
from evolution import *
import numpy as np
import matplotlib.pyplot as plt

TT = TexTab()
EM = EvolutionaryModels(load_tr_pickle = True)

#### DATA
inputf = ascii.read(args.filename, header_start=0, delimiter=",")
inputf['RA'], inputf['DEC'] = to_deg(inputf['RA'],inputf['DEC'])
input_tab = Table({'STAR' : inputf['STAR'], 'RA' : inputf['RA'], 'DEC' : inputf['DEC'], 'SBTYPE' : inputf['SBTYPE'], 'SSPOC': inputf['SSPOC'], 'SFFI' : inputf['SFFI'], 'MASK' : inputf['MASK'], 'SEDFIT' : inputf['SEDFIT'], 'RDMAG' : inputf['RDMAG'], 'VART' : inputf['VART']})

###### CROSS-MATCHING AND DATA PREPARATION
data = XMatching(input_tab, args.xmtabs, args.xmcols, vizier = args.xmviz, load_pickle = True).matched
data['EDD'] = (KAPPA * SBOLTZ / Ccm) * (data['TEFF']**4) / 10**data['LOGG']
data['MDOT'] = np.log10(1e-6*data['MDOT'])
data['TEFF'] = np.ma.log10(data['TEFF'])
#TT.TabSample(data)

#print data['STAR','RA','DEC','TIC','SSPOC','SFFI'].pprint(max_lines = -1, max_width = -1)
'''data['MEVOL'] = EM.interp2d('logTe', 'logg', 'Mass', data['TEFF'], data['LOGG'], post_rsg = False, method='linear')
EM.plot_spectroHR(data, output_format = 'eps', output_name = 's_hrdiag', hold = True, inter = False)
x_range = np.arange(3.6,4.9,0.1)
f = lambda x, gamma : 4 * x - np.ma.log10(gamma / (KAPPA * SBOLTZ / Ccm))
g1, = EM.ax.plot(x_range,f(x_range, 1),'r-', lw = 2)
g06, = EM.ax.plot(x_range,f(x_range, 0.6),'r:', lw = 2)
EM.ax.legend([g1, g06], [ r'$\Gamma = 1$',r'$\Gamma = 0.6$'])
EM.panel.PanelClose()'''


###### EXTRACTION

# LIGHTCURVES FROM TPF - BREAK INTO 3 OR 4 PARTS DEPENDING ON THE TYPES
#ExtractLC(data[:20], 5, all_type = False, output_format = None, save_files = False, inter=False, output_name = 'XLC_v2')

# PERIODOGRAMS AND RED NOISE MODELS - PREWHITENING - ORGANIZING FREQUENCIES - BREAK INTO 3 PARTS
# For published prewhitening figure, plot data[24:25] in 'eps'
#ExtractLS(data,  prewhiten = True, output_format = None, output_name = 'XLS_v2', rows_page = 6, cols_page = 2, inter = False)
#import os; os.system('systemctl poweroff')

FM = FrequencyManager(data)


data['LOGG'][data['STAR'] == 'HD152236'] = 1.5
###### COLLECT PHOTOMETRIC AND SYNTHETIC DATA
photo = PhotoCollector(data, load_pickle = True)
photo_data = photo.photo_tab 
filter_dict = photo.filt_dict
model = model_dict(filter_dict, load_pickle = True)
#TT.TabPhoto(photo_data)

###### SET FITTING AND VISUALIZING - BUILD/LOAD TABLE WITH SED PARAMETERS
#SEDBuilder(photo_data[:14], filter_dict, fit_sed = True, fit_model_dict = model, fit_bodies = 'p',
#		 rows_page = 7, cols_page = 2, output_name = 'SED', coll_x = True, coll_y = True, output_format = 'eps', 
#		 figsize = (10,18))
#SEDBuilder(photo_data[14:], filter_dict, fit_sed = True, fit_model_dict = model, fit_bodies = 'p',
#		 rows_page = 8, cols_page = 3, output_name = 'SED', coll_x = True, coll_y = True, output_format = 'eps', 
#		 figsize = (13,16))

SED = SEDBuilder(photo_data, filter_dict, fit_sed = True, fit_model_dict = model, fit_bodies = 'p',
		output_format = None, load_pickle = True)
data['S_MASS'] = np.ma.log10(radmass(data['LOGG'],SED.sed_tab['S_RAD']))
data['GABS'] = photo_data['Gmag'] - 5*np.ma.log10(SED.sed_tab['S_DIST']) + 5 - SED.sed_tab['A_V'] 
data['S_RAD'] = SED.sed_tab['S_RAD']
#TT.TabSED(SED.sed_tab)

data['VCRIT'] = EM.interp2d('logTe', 'logL', 'vcrit1', data['TEFF'], SED.sed_tab['S_LOGL'],method='linear')
data['VCRIT'].format = '.1f'  #;  EM.plot('logTe', 'logL', 'vcrit1', levels = 30)

###### PROCESSING/VISUALIZING LIGHTCURVES/PERIODOGRAMS

# LC - CALCULATE/LOAD METRICS
LC = Processing(data, level='lc', output_format = None, load_mt_pickle = True)

# RN - CALCULATE/LOAD METRICS
LS = Processing(hstack([data,FM.freq_tab['FF','A_FF']]), level='ls', output_format = None, load_rn_pickle = True)
LS.rn_tab['VCHAR'] = np.log10( 1 / (2 * PI * LS.rn_tab['TAU']) )
data['F_ROT'] = LS.rot_stat
data['INDFF'] = LS.indep_freq
data['INDFFS'] = LS.indep_freq_sec

 
# FOR VISUALIZING
#Processing(data[:15], level='lc', rows_page = 5, cols_page = 3, figsize = (16,10),
#			  output_name = 'LC', output_format = 'eps', inter = False)
#Processing(data[15:], level='lc', rows_page = 8, cols_page = 3, figsize = (18,18),
#			  output_name = 'LC', output_format = 'eps', inter = False)
#Processing(data[:14], level='ls', rows_page = 7, cols_page = 2, figsize = (10,18),
#			  output_name = 'LS', output_format = 'eps', inter = False)
#Processing(data[14:], level='ls', rows_page = 8, cols_page = 3, figsize = (13,16),
#			  output_name = 'LS', output_format = 'eps', inter = False)


#TT.TabCalcProp(hstack([SED.sed_tab['STAR','A_V','LOGC','S_LOGL'], LC.mt_tab['SVAR','ETA','PSI','SKEW'],
#	LS.rn_tab['LOGW','LOGR0','TAU','GAMMA']]))


###### STATISTICS
corr_tab = hstack([data,LC.mt_tab,SED.sed_tab,LS.rn_tab,FM.freq_tab])
print corr_tab['STAR_1','SVAR','ETA','S_LOGL','SKEW','PSI','INDFF','INDFFS','VART','F_ROT'][corr_tab['F_ROT']==False].pprint(max_lines=-1)

# CORRELATION BETWEEN STELLAR PARAMETERS AND METRICS
#corr_scatter(corr_tab,x = ['TEFF','S_LOGL','EDD'], y = ['SVAR','ETA','PSI','SKEW'],
#	     mode = 'matrix', output_name = 'scatter_par', output_format = 'eps', figsize = (11,10), inter = False)
#corr_scatter(corr_tab,x = ['NABUN','VMIC','VMAC','VINF','BETA'], y = ['SVAR','ETA','PSI','SKEW'],
#	     mode = 'matrix', output_name = 'scatter_par', output_format = 'eps', figsize = (11,12), inter = False)

# CORRELATION BETWEEN STELLAR AND RED NOISE PARAMETERS
'''
corr_scatter(corr_tab, x = ['TEFF','S_LOGL','EDD'],y = ['LOGW','LOGR0','TAU','GAMMA'], 
	    figsize = (11,10), output_name = 'scatter_rn', output_format = 'eps', inter = False)
corr_scatter(corr_tab, x = ['NABUN','VMIC','VMAC','VINF','BETA'],y = ['LOGW','LOGR0','TAU','GAMMA'], 
	    figsize = (11,12), output_name = 'scatter_rn', output_format = 'eps', inter = False)



# COMPARISON TO BOWMAN+19
cr = corr_scatter(corr_tab,x = ['zero','VCHAR'], y = ['GABS'],
	     mode = 'matrix', rows_page = 2, cols_page = 1, output_format = None, output_name = 'rn_gmag', figsize = (7,11), hold=True, inter = False)
p1 = cr.get_ax(0,0)
p2 = cr.get_ax(1,0)
xr = np.array([-2,-12])
p1.invert_xaxis(); del p1.texts[0]; 
B19A = lambda x: -0.2 * x + 0.3 - 3 - 2.2
b19A_line, = p1.plot(xr,B19A(xr),'b--')
p1.set_ylim(-4.99,-2.1)
p1.legend([b19A_line], [r'$-$0.2x$-$2.7 ($-$const.) (Bowman et al. 2019)'])#,loc=3)
p2.invert_xaxis(); del p2.texts[0]
mask_ind = (corr_tab['GABS'].mask == False) & (corr_tab['VCHAR'].mask == False)
par2 = np.polyfit(corr_tab['GABS'][mask_ind], corr_tab['VCHAR'][mask_ind], 1)
M2 = np.poly1d(par2)
B19B = lambda x: 0.16 * x + 0.06 + 1.2  # Added C2 = 1.2
b19B_line, = p2.plot(xr,B19B(xr),'b--')
fit_line, = p2.plot(xr,M2(xr),'r',lw=1.5)
p2.set_ylim(-0.39,0.9)
p2.legend([fit_line, b19B_line], [r'%.2fx+%.2f' % tuple(par2), r'0.16x+0.06 (+const.) (Bowman et al. 2019)'])#,loc=3)
cr.GridClose()
'''
'''
# CORRELATION BETWEEN STELLAR PARAMETERS AND FREQUENCIES
#corr_scatter(corr_tab,x = ['TEFF','S_LOGL'], y = ['FF','A_FF'],
#	     mode = 'frequency', output_format = 'eps', output_name = 'freq_scatter', figsize = (9,11), inter = False)

FH = freq_hist(corr_tab, x = ['TEFF','S_LOGL'], y = ['FF'], snr_thres = 0,  figsize = (9,11), output_format = None, output_name = 'freq_hist', ngroups = 3, inter = False)


###### STELLAR EVOLUTION
HR = HRdiagram(hstack([data,SED.sed_tab,FH.outl_tab]), tkey = 'TEFF', lkey = 'S_LOGL', figsize = (9,12), output_format = 'eps', output_name = 'hrdiag', hold = True, inter = False)
#for i,j in zip(FH.thres_group['TEFF'],FH.thres_group['S_LOGL']) :
#		HR.ax.axhline(y=j, color='c', ls='-')
#		HR.ax.axvline(x=i, color='c', ls='-')
HR.ax.plot([3.8,4.9],[np.log10((18**2)*((10**x)/5778)**4) for x in [3.8,4.9]],'c:')
HR.ax.text(4.12,3.98,r'18 R$_\odot$', weight= 'bold')
HR.ax.plot([3.8,4.6],[np.log10((50**2)*((10**x)/5778)**4) for x in [3.8,4.6]],'m:')
HR.ax.text(3.98,4.32,r'50 R$_\odot$', weight= 'bold')
HR.ax.set_xlim(4.77,3.9)
HR.ax.set_ylim(3.9,6.4)
HR.PanelClose()'''
