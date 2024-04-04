#%%
import args
from astropy.table import Table
from astroquery.xmatch import XMatch
from astropy.io import ascii
from functions import *
from time_series import *
from statistics import *
from sed import *
#from tables import *
from constants import *
from evolution import *
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as plt

TT = TexTab()

#### DATA
inputf = ascii.read(args.filename, header_start=0, delimiter=",")
inputf['RA'], inputf['DEC'] = to_deg(inputf['RA'],inputf['DEC'])
input_tab = Table({'STAR' : inputf['STAR'], 'RA' : inputf['RA'], 'DEC' : inputf['DEC'], 'SBTYPE' : inputf['SBTYPE'], 'SSPOC': inputf['SSPOC'], 'SFFI' : inputf['SFFI'], 'MASK' : inputf['MASK'], 'SEDFIT' : inputf['SEDFIT'], 'RDMAG' : inputf['RDMAG']})

# CROSS-MATCHING AND DATA PREPARATION
data = XMatching(input_tab, args.xmtabs, args.xmcols, vizier = args.xmviz, load_pickle = True).matched
data['EDD'] = (KAPPA * SBOLTZ / Ccm) * (data['TEFF']**4) / 10**data['LOGG']
data['MDOT'] = np.log10(1e-6*data['MDOT'])
data['TEFF'] = np.log10(data['TEFF'])
#TT.TabSample(data)


# EXTRACTION OF TESS LC AND LS
# For ALL_TYPES BREAK INTO 8 (BRIGHT) - 30 - 55 - end
# For SELECTED TYPES BREAK INTO 8 (BRIGHT) - 38 - end
#### EXTRACT LIGHTCURVES
#ExtractLC(data, 4, all_type = False, output_format = 'pdf', save_files = False, inter=False, output_name = 'XLC_v2')

#### EXTRACT PERIODOGRAMS AND RED NOISE MODELS - PERFORM PREWHITENING
# For published prewhitening figure, plot data[24:25] in 'eps'
# BREAK INTO 20 - 40 - end
#ExtractLS(data,  prewhiten = True, output_format = None, output_name = 'XLS_v2', rows_page = 6, cols_page = 2, inter = False)

#os.system('systemctl poweroff') 

#### Processing LIGHTCURVES OR LOMB-SCARGLE - 
#Processing(data[:15], level='lc', rows_page = 5, cols_page = 3, figsize = (16,10),
#			  output_name = 'LC', output_format = 'eps', inter = False)
#Processing(data[15:], level='lc', rows_page = 8, cols_page = 3, figsize = (18,18),
#			  output_name = 'LC', output_format = 'eps', inter = False)
#Processing(data[:14], level='ls', rows_page = 7, cols_page = 2, figsize = (10,18),
#			  output_name = 'LS', output_format = 'eps', inter = False)
#Processing(data[14:], level='ls', rows_page = 8, cols_page = 3, figsize = (13,16),
#			  output_name = 'LS', output_format = 'eps', inter = False)

#### SAVE OR LOAD LC METRICS TABLE
LC = Processing(data, level='lc', output_format = None, load_mt_pickle = True)
LC.mt_tab['e_ETA'] = LC.mt_tab['e_ETA'] /(np.log(10) * (LC.mt_tab['ETA'])); LC.mt_tab['ETA'] = np.log10(LC.mt_tab['ETA']); 


#### SAVE OR LOAD RED NOISE PARAMETER TABLE
LS = Processing(data, level='ls', output_format = None, load_rn_pickle = True)
LS.rn_tab['VCHAR'] = np.log10( 1 / (2 * PI * LS.rn_tab['tau']) )
#TT.TabRN(LS.rn_tab)
LS.rn_tab['tau'] = np.log10(LS.rn_tab['tau'])


###### PHOTOMETRIC AND SYNTHETIC DATA
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
SED.sed_tab['S_MASS'] = radmass(data['LOGG'],SED.sed_tab['S_RAD'])
SED.sed_tab['GABS'] = photo_data['Gmag'] - 5*np.log10(SED.sed_tab['S_DIST']) + 5 - SED.sed_tab['A_V'] 
#TT.TabSED(SED.sed_tab)


###### STATISTICS
corr_tab = hstack([data,LC.mt_tab,SED.sed_tab,LS.rn_tab])
corr_tab['EDD'][corr_tab['STAR_1'] == 'HD152236'] = np.nan

# 1 - Correlations between stellar parameters and metrics
corr_scatter(data = corr_tab,x = ['MASS','S_LOGL','EDD'], y = ['SVAR','ZCROSS','PSI','ETA'],
	     mode = 'matrix', rows_page = 2, cols_page = 2, output_format = None, figsize = (15,8), inter = True)
'''corr_scatter(data = corr_tab,x = ['TEFF','NABUN','VMIC','VMAC'], y = ['SVAR','ZCROSS','PSI','ETA'],
	     mode = 'matrix', output_format = None, figsize = (15,8), inter = True)
corr_scatter(data = corr_tab,x = ['MASS','S_LOGL','EDD'], y = ['w','zero','tau','gamma'],
	     mode = 'matrix', output_format = None, figsize = (15,8), inter = True)
corr_scatter(data = corr_tab,x = ['TEFF','NABUN','VMIC','VMAC'], y = ['w','zero','tau','gamma'],
	     mode = 'matrix', ouutput_format = None, figsize = (15,8), inter = True)'''

# 2 - Correlations between stellar parameters and red noise parameters
'''cr = corr_scatter(jtab = corr_tab,
	     x = ['MASS','S_LOGL','EDD','MDOT'],
	     #x = ['TEFF','NABUN','VMIC','VMAC'],
	     y = ['SVAR','ZCROSS','PSI','ETA'],
	     coll_x = True, coll_y = True, output_format = None, hold=True, inter = True)'''

# 3 - Comparison to Bowman+19
cr = corr_scatter(data = corr_tab,x = ['zero','VCHAR'], y = ['GABS'],
	     mode = 'matrix', rows_page = 2, cols_page = 1, output_format = 'eps', output_name = 'rn_gmag', figsize = (7,11), hold=True, inter = False)
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


# 4 - Correlations between stellar parameters and distribution of frequencies
#TH = tess_hist(ext_tab = photo_data[ext_keys], snr_show = 4, output_format = None, coll_x = True, inter = True)
#photo_data = hstack([photo_data,TH.hist_tab])
#HRdiagram(photo_data, lkey = 'LUM', tlines = TH.prop_thres['TEFF'], llines = TH.prop_thres['LUM'], inter = True)


# 5 - Histograms





#HRdiagram(photo_data, cbar_key = 'MASS', lkey = 'LUM', inter = True)










