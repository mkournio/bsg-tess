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

TT = TexTab()

#### DATA
inputf = ascii.read(args.filename, header_start=0, delimiter=",")
inputf['RA'], inputf['DEC'] = to_deg(inputf['RA'],inputf['DEC'])
input_tab = Table({'STAR' : inputf['STAR'], 'RA' : inputf['RA'], 'DEC' : inputf['DEC'], 'SBTYPE' : inputf['SBTYPE'], 'SSPOC': inputf['SSPOC'], 'SFFI' : inputf['SFFI'], 'MASK' : inputf['MASK'], 'SEDFIT' : inputf['SEDFIT'], 'RDMAG' : inputf['RDMAG']})

# CROSS-MATCHING
data = XMatching(input_tab, args.xmtabs, args.xmcols, vizier = args.xmviz, load_pickle = True).matched
#TT.TabSample(data)


# For ALL_TYPES BREAK INTO 8 (BRIGHT) - 30 - 55 - end
# For SELECTED TYPES BREAK INTO 8 (BRIGHT) - 38 - end
#### EXTRACT LIGHTCURVES
#ExtractLC(data, 4, all_type = False, output_format = 'pdf', save_files = False, inter=False, output_name = 'XLC_v2')

#### EXTRACT PERIODOGRAMS AND RED NOISE MODELS - PERFORM PREWHITENING
# For published prewhitening figure, plot data[24:25] in 'eps'
# BREAK INTO 20 - 40 - end
#ExtractLS(data,  prewhiten = True, output_format = None, output_name = 'XLS_v2', rows_page = 6, cols_page = 2, inter = False)

#os.system('systemctl poweroff') 

#### VISUALIZE LIGHTCURVES OR LOMB-SCARGLE - 
#Visualize(data[:15], level='lc', rows_page = 5, cols_page = 3, figsize = (16,10),
#			  output_name = 'LC', output_format = 'eps', inter = False)
#Visualize(data[15:], level='lc', rows_page = 8, cols_page = 3, figsize = (18,18),
#			  output_name = 'LC', output_format = 'eps', inter = False)
#Visualize(data[:15], level='ls', rows_page = 5, cols_page = 3, figsize = (15,12),
#			  output_name = 'LS', output_format = 'eps', inter = False)
#Visualize(data[15:], level='ls', rows_page = 8, cols_page = 3, figsize = (18,20),
#			  output_name = 'LS', output_format = 'eps', inter = False)

#### LOAD RED NOISE PARAMETER TABLE
LS = Visualize(data, level='ls', output_format = None, load_rn_pickle = True)


###### MATCH DATA TABLE WITH PHOTOMETRY.
photo = PhotoCollector(data, load_pickle = True)
photo_data = photo.photo_tab 
filter_dict = photo.filt_dict

#TT.TabPhoto(photo_data)


###### BUILD/LOAD TABLE WITH SYNTHETIC SED FLUXES
model = model_dict(filter_dict, load_pickle = True)

###### FIT SEDS AND VISUALIZE - BUILD/LOAD TABLE WITH SED PARAMETERS
SED = SEDBuilder(photo_data, filter_dict, fit_sed = True, fit_model_dict = model, fit_bodies = 'p',
		 rows_page = 4, cols_page = 5, output_name = 'SED_v1', coll_x = True, output_format = 'pdf', 
		 inter = False, load_pickle = False)

TT.TabSED(hstack([LS.rn_tab,SED.sed_tab]))


#HRdiagram(sed, lkey = 'LUM', output_format = None, inter = True)

#%%


# STATISTICS
#ext_keys = ['STAR','TEFF','NABUN','LOGL','LUM','MASS','VMIC']
# Correlation matrix for the studied parameters to enable feature selection
#autocorr_scatter(vector=sed[ext_keys], output_format = None, inter = True, coll_x = True, coll_y = True )


# Correlations between stellar parameters and red noise
#corr_scatter(ext_x = rn, ext_y = sed[ext_keys], match_keys = 'STAR', 
#             coll_x = True, coll_y = True, output_format = 'pdf', inter = False)


#Correlations between stellar parameters and distribution of frequencies
#TH = tess_hist(ext_tab = photo_data[ext_keys], snr_show = 4, output_format = None, coll_x = True, inter = True)
#photo_data = hstack([photo_data,TH.hist_tab])
#HRdiagram(photo_data, lkey = 'LUM', tlines = TH.prop_thres['TEFF'], llines = TH.prop_thres['LUM'], inter = True)


#HRdiagram(photo_data, cbar_key = 'MASS', lkey = 'LUM', inter = True)


######## NOTES

#ext_keys = ['TEFF','VSINI','MDOT','LOGLM','LOGQ']
#ext_keys = ['TEFF','VSINI','MDOT','LOGLM','LOGQ','NABUN','LOGD','LOGG','MASS','LOGL','VMAC','VMIC']
#ext_keys = ['TEFF','MDOT','VSINI','LOGQ','NABUN','LOGD','LOGLM']
#ext_keys = ['LOGG','MASS','LOGL','VMAC','VMIC']
#ext_keys = ['EW3995','EW4129','EW4131','EW4553','EW4568','EW4575']








