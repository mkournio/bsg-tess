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
try:
	data = pickle.load(open(PICKLE_PATH+'input.pkl','rb'))
	print 'Loaded pickle: input'

except (OSError, IOError) as e:
	# READING INPUT
	inputf = ascii.read(args.filename, header_start=0, delimiter=",", data_start=args.from_line, data_end=args.to_line)
	inputf['RA'], inputf['DEC'] = to_deg(inputf['RA'],inputf['DEC'])
	input_tab = Table({'STAR' : inputf['STAR'], 'RA' : inputf['RA'], 'DEC' : inputf['DEC'], 'TMAG' :  inputf['TMAG'] , 'SBTYPE' : inputf['SBTYPE'], 'SSPOC': inputf['SSPOC'], 'SFFI' : inputf['SFFI'], 'MASK' : inputf['MASK'], 'SEDFIT' : inputf['SEDFIT'], 'RDMAG' : inputf['RDMAG'],})

	# CROSS-MATCHING
	xm = XMatching(input_tab)
	xm.match(xtabs = ['HAUCKE+19','FRASER+10'], xcols = ['TEFF','NABUN','MASS','LOGL','VMIC', 'LOGG'])
	xm.match(xtabs = ['vizier:I/347/gaia2dis'], xcols = ['DIST'], viz_col = 'rest')
	xm.match(xtabs = ['vizier:IV/39/tic82'], xcols = ['TIC'], viz_col = 'TIC', column_format = 'int')
	data = xm.matched
	pickle.dump(data,open(PICKLE_PATH+'input.pkl','wb'))
	print 'Task executed: input'

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

#### VISUALIZE LIGHTCURVES OR LOMB-SCARGLE - BUILD RED NOISE PARAMETER TABLE
# For publishing, break into 24 - 48 - end
#Visualize(data[48:], level='lc', rows_page = 8, cols_page = 3, 
#			  output_name = 'LC', output_format = 'eps', inter = False)

#Visualize(data[:10], level='ls', rows_page = 5, cols_page = 2, 
#			  output_name = 'LS', output_format = 'eps', inter = False)




try:
	rn = pickle.load(open(PICKLE_PATH+'rn.pkl','rb'))
	print 'Loaded pickle: red noise properties'
except:	
	# Visualize light curves OR power spectra
	LS = Visualize(data, level='ls', rows_page = 5, cols_page = 2, 
			  output_name = 'LS', output_format = None, inter = False)
	rn = LS.rn_tab
	pickle.dump(rn,open(PICKLE_PATH+'rn.pkl','wb'))	
	print 'Task executed: input'

print rn
'''

###### MATCH DATA TABLE WITH PHOTOMETRY
try:
	photo = pickle.load(open(PICKLE_PATH+'photo.pkl','rb'))
except (OSError, IOError) as e:
	# COLLECT PHOTOMETRY FROM LITERATURE
	photo = PhotoCollector(data)
	pickle.dump(photo,open(PICKLE_PATH+'photo.pkl','wb'))

photo_data = photo.photo_tab #; photo.save_tab()
filter_dict = photo.filt_dict


###### BUILD TABLE WITH SYNTHETIC SED FLUXES
try:
	model = pickle.load(open(PICKLE_PATH+'model.pkl','rb'))	
except (OSError, IOError) as e:
	# CREATE SED MODEL TABS (Kurucz & PoWR)
	synth_l = []
	for v in filter_dict.values(): synth_l += v['l']; synth_l = sorted(set(synth_l))
	KuruczTab = synthetic_fluxes(synth_l,model_type = 'kurucz')
	PoWRTab = synthetic_fluxes(synth_l,model_type = 'powr')
	model = {'powr' : PoWRTab, 'kurucz' : KuruczTab}
	pickle.dump(model,open(PICKLE_PATH+'model.pkl','wb'))


###### FIT SEDS AND VISUALIZE - BUILD TABLE WITH SED PARAMETERS
try:
	sed = pickle.load(open(PICKLE_PATH+'sed.pkl','rb'))	
except (OSError, IOError) as e:
	# BUILD SPECTRAL ENERGY DISTRIBUTION - APPEND TO TABLE
	SED = SEDBuilder(photo_data, filter_dict, fit_sed = True, fit_model_dict = model,
		 rows_page = 4, cols_page = 5, output_name = args.filename+'_SED', 
		 coll_x = True, output_format = None, inter = False)
	sed = hstack([photo_data,SED.fit_tab])
	pickle.dump(sed,open(PICKLE_PATH+'sed.pkl','wb'))

HRdiagram(sed, lkey = 'LUM', output_format = None, inter = True)

#%%


# STATISTICS
ext_keys = ['STAR','TEFF','NABUN','LOGL','LUM','MASS','VMIC']
# Correlation matrix for the studied parameters to enable feature selection
#autocorr_scatter(vector=sed[ext_keys], output_format = None, inter = True, coll_x = True, coll_y = True )


# Correlations between stellar parameters and red noise
corr_scatter(ext_x = rn, ext_y = sed[ext_keys], match_keys = 'STAR', 
             coll_x = True, coll_y = True, output_format = 'pdf', inter = False)


#Correlations between stellar parameters and distribution of frequencies
TH = tess_hist(ext_tab = photo_data[ext_keys], snr_show = 4, output_format = None, coll_x = True, inter = True)
photo_data = hstack([photo_data,TH.hist_tab])
HRdiagram(photo_data, lkey = 'LUM', tlines = TH.prop_thres['TEFF'], llines = TH.prop_thres['LUM'], inter = True)


HRdiagram(photo_data, cbar_key = 'MASS', lkey = 'LUM', inter = True)


######## NOTES

#ext_keys = ['TEFF','VSINI','MDOT','LOGLM','LOGQ']
#ext_keys = ['TEFF','VSINI','MDOT','LOGLM','LOGQ','NABUN','LOGD','LOGG','MASS','LOGL','VMAC','VMIC']
#ext_keys = ['TEFF','MDOT','VSINI','LOGQ','NABUN','LOGD','LOGLM']
#ext_keys = ['LOGG','MASS','LOGL','VMAC','VMIC']
#ext_keys = ['EW3995','EW4129','EW4131','EW4553','EW4568','EW4575']


'''






