import argparse
from astropy.table import Table
from astroquery.xmatch import XMatch
from astropy.io import ascii
from functions import *
from time_series import *
from statistics import *

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename",default= None, help="File containing the star table")
parser.add_argument("-l1","--from_line", type=int, default = 1, help="First reading line")
parser.add_argument("-l2","--to_line", type=int, default = None, help="End reading line")	
args = parser.parse_args()
	
	#### READ TABLE
inputf = ascii.read(args.filename, header_start=0, delimiter=",",data_start=args.from_line, data_end=args.to_line)
inputf['RA'], inputf['DEC'] = to_deg(inputf['RA'],inputf['DEC'])

coord = Table({'STAR' : inputf['STAR'], 'RA' : inputf['RA'], 'DEC' : inputf['DEC'], 
		       'SSPOC': inputf['SSPOC'], 'SFFI' : inputf['SFFI'], 'MASK' : inputf['MASK']})


data = XMatch.query(cat1=coord, cat2='vizier:IV/39/tic82',  max_distance=2*u.arcsec, colRA1='RA', colDec1='DEC')
fdata =  match_tabs(coord,data)
fdata.pprint(max_lines=-1)
	
	#### EXTRACT LIGHTCURVES (ALL)
	#ExtLCs(data['STAR'],data['RA'],data['DEC'],data['TIC'],4,pdf_name=args.filename+'_ExtLCs')

	#### EXTRACT LIGHTCURVES (SELECTED)
	#ExtLCs(data['STAR'],data['RA'], data['DEC'], data['TIC'], 4, sspoc = data['SSPOC'], sffi = data['SFFI'],
	#		    thmask = data['MASK'], pdf_name=args.filename+'_ExtLCs')

	#### VISUALIZE LIGHTCURVES
	#VisLCs(star_ids=fdata['STAR'], plots_per_grid = 8, folder = 'EXTLC_R1', pdf_name = args.filename+'_LCs', 
	#								       pdf_save = True, inter = False)

	#### VISUALIZE PERIODOGRAMS
	#ProcessLevel(level='ls',star_ids=fdata['STAR'], rows_page = 6, cols_page = 2, folder = 'EXTLC_R1',
	#		        output_name = args.filename+'_LC', output_format = None, inter = False)

LS = Visualize(level='ls',star_ids=fdata['STAR'], rows_page = 5, cols_page = 5, 
			   output_name = args.filename+'_LS', coll_x = True, coll_y = True, 
			output_format = None, inter = False)

ext_files = ['HAUCKE+19','FRASER+10']
ext_keys = ['TEFF','MDOT','VSINI','LOGQ','NABUN','LOGD','LOGLM']
#ext_keys = ['LOGG','MASS','LOGL','VMAC','VMIC','VINF']
#ext_keys = ['EW3995','EW4129','EW4131','EW4553','EW4568','EW4575']
for k in ext_keys:
	external = XmExtCol(inputf['RA'], inputf['DEC'],ext_files=ext_files,ext_col= k)
	fdata = hstack([fdata,external])

ext_keys = ['STAR'] + ['f_'+ k for k in ext_keys] + ext_keys


Correlation(ext_x = LS.rn_tab, ext_y = fdata[ext_keys], match_keys = 'STAR', 
		coll_x = True, coll_y = True, output_format = None, inter = True)



















