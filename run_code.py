import argparse
from astropy.table import Table
from astroquery.xmatch import XMatch
from astropy.io import ascii
from functions import *
from time_series import *
from statistics import *
from sed import *
#from tables import *

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename",default= None, help="File containing the star table")
parser.add_argument("-l1","--from_line", type=int, default = 1, help="First reading line")
parser.add_argument("-l2","--to_line", type=int, default = None, help="End reading line")	
args = parser.parse_args()
	
# READING INPUT
inputf = ascii.read(args.filename, header_start=0, delimiter=",",data_start=args.from_line, data_end=args.to_line)
inputf['RA'], inputf['DEC'] = to_deg(inputf['RA'],inputf['DEC'])
input_tab = Table({'STAR' : inputf['STAR'], 'RA' : inputf['RA'], 'DEC' : inputf['DEC'], 
		       'SSPOC': inputf['SSPOC'], 'SFFI' : inputf['SFFI'], 'MASK' : inputf['MASK'], 'SEDFIT' : inputf['SEDFIT']})


# CROSS-MATCHING
xm = XMatching(input_tab)
xm.match(xtabs = ['HAUCKE+19','FRASER+10'], xcols = ['TEFF','NABUN','MASS','LOGL','VMIC', 'LOGG'])
xm.match(xtabs = ['vizier:I/347/gaia2dis'], xcols = ['DIST'], viz_col = 'rest')
xm.match(xtabs = ['vizier:IV/39/tic82'], xcols = ['TIC'], viz_col = 'TIC', column_format = 'int')
data = xm.matched
#data.pprint(max_lines=-1)


'''
#### EXTRACT LIGHTCURVES (ALL)
#ExtLCs(data['STAR'],data['RA'],data['DEC'],data['TIC'],4,pdf_name=args.filename+'_ExtLCs')

#### EXTRACT LIGHTCURVES (SELECTED)
#ExtLCs(data['STAR'],data['RA'], data['DEC'], data['TIC'], 4, sspoc = data['SSPOC'], sffi = data['SFFI'],
#		    thmask = data['MASK'], pdf_name=args.filename+'_ExtLCs')


#### VISUALIZE PERIODOGRAMS
#ProcessLevel(level='ls',star_ids=fdata['STAR'], rows_page = 6, cols_page = 2, folder = 'EXTLC_R1',
#		        output_name = args.filename+'_LC', output_format = None, inter = False)


# Visualize light curves OR power spectra
#LS = Visualize(level='ls',star_ids=fdata['STAR'], rows_page = 5, cols_page = 5, 
#			  output_name = args.filename+'_LS', coll_x = True, coll_y = True, 
#			  output_format = None, inter = False)
'''

# COLLECT PHOTOMETRY FROM LITERATURE
photo = PhotoCollector(data)
photo_data = photo.photo_tab#; photo.save_tab()
filter_dict = photo.filt_dict


# CREATE SED MODEL TABS (Kurucz & PoWR)
synth_l = []
for v in filter_dict.values(): synth_l += v['l']; synth_l = sorted(set(synth_l))
KuruczTab = synthetic_fluxes(synth_l,model_type = 'kurucz')
PoWRTab = synthetic_fluxes(synth_l,model_type = 'powr')
model_dict = {'powr' : PoWRTab, 'kurucz' : KuruczTab}

# BUILD SPECTRAL ENERGY DISTRIBUTION
SED = SEDBuilder(photo_data, filter_dict, fit_sed = True, fit_dict = {'MODELS': model_dict, 'FIT_BODIES': ['psph']},
		 rows_page = 4, cols_page = 5, output_name = args.filename+'_SED', 
		 coll_x = True, output_format = None, inter = True)

'''
# Correlation matrix for the studied parameters to enable feature selection
#autocorr_scatter(vector=fdata[ext_keys], output_format = None, inter = 'True', coll_x = True, coll_y = True )

# Correlations between stellar parameters and red noise
#corr_scatter(ext_x = LS.rn_tab, ext_y = fdata[ext_keys], match_keys = 'STAR', 
#             coll_x = True, coll_y = True, output_format = 'pdf', inter = False)

#Correlations between stellar parameters and distribution of frequencies

tess_hist(ext_tab = fdata[ext_keys], snr_show = 4, output_format = None, coll_x = True, inter = True)

'''






######## NOTES

#ext_keys = ['TEFF','VSINI','MDOT','LOGLM','LOGQ']
#ext_keys = ['TEFF','VSINI','MDOT','LOGLM','LOGQ','NABUN','LOGD','LOGG','MASS','LOGL','VMAC','VMIC']
#ext_keys = ['TEFF','MDOT','VSINI','LOGQ','NABUN','LOGD','LOGLM']
#ext_keys = ['LOGG','MASS','LOGL','VMAC','VMIC']
#ext_keys = ['EW3995','EW4129','EW4131','EW4553','EW4568','EW4575']









