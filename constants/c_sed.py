#### Constants for spectral energy distributions

#Photometric catalogs
PHOTO_CATS =	{
	'Merm91':	{# Homogeneous Means in the UBV System (Mermilliod 1991)
			'act'	: 1,		 
			'fit'	: 1,
			'cat'	: 'II/168/ubvmeans', 
			'rad'	: 4, 
			'mrk'	: {'m': 'x', 'c':'b', 's':9}, 
			'flt'	: ['Vmag','B-V','U-B','e_Vmag','e_B-V','e_U-B']
			}, 
	'NOMAD05':	{# NOMAD Catalog (Zacharias+ 2005)
			'act'	: 1,		
			'fit'	: 1,
			'cat'	: 'I/297/out',
			'rad'	: 2,
			'mrk'	: {'m': '+', 'c':'b', 's':9}, 
			'flt'	: ['Bmag','Vmag','Rmag']
			},
	'APASS':	{# AAVSO Photometric All Sky Survey (APASS) DR9 (Henden+, 2016)
			'act'	: 1,		
			'fit'	: 1,
			'cat'	: 'II/336/apass9',
			'rad'	: 1,
			'mrk'	: {'m': '.', 'c':'b', 's':9}, 
			'flt'	: ['Bmag','Vmag','e_Bmag','e_Vmag','gpmag','e_gpmag','rpmag','e_rpmag','ipmag','e_ipmag']
			},
	'Pan-STARRS':	{# The Pan-STARRS release 1 (PS1) Survey - DR1 (Chambers+, 2016)	
			'act'	: 0,		 
			'fit'	: 0,
			'cat'	: 'II/349/ps1',
			'rad'	: 5,
			'mrk'	: {'m': 'v', 'c':'c', 's':9}, 
			'flt'	: ['gmag','rmag','imag','zmag','ymag','e_gmag','e_rmag','e_imag','e_zmag','e_ymag']
			},
	'GAIA2':	{# Gaia DR2 (Gaia Collaboration, 2018)
			'act'	: 0,				
			'fit'	: 0,
			'cat'	: 'I/345/gaia2',
			'rad'	: 2,
			'mrk'	: {'m': 'o', 'c':'r', 's':9},
			'flt'	: ['phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag']
			},
	'GAIA3':	{# Gaia DR3 (Gaia Collaboration, 2022)	
			'act'	: 1,				
			'fit'	: 1,
			'cat'	: 'I/355/gaiadr3',
			'rad'	: 2,
			'mrk'	: {'m': 'o', 'c':'k', 's':9},
			'flt'	: ['BPmag','Gmag','RPmag','e_BPmag','e_Gmag','e_RPmag']
			},
	'2MASS':	{# 2MASS All-Sky Catalog of Point Sources (Cutri+ 2003)	
			'act'	: 1,		
			'fit'	: 1,
			'cat'	: 'II/246/out',
			'rad'	: 2,
			'mrk'	: {'m': '^', 'c':'g', 's':9},
			'flt'	: ['Jmag','Hmag','Kmag','e_Jmag','e_Hmag','e_Kmag']
			},
	'GLIMPSE':	{# GLIMPSE Source Catalog (I + II + 3D) (IPAC 2008)
			'act'	: 0,		
			'fit'	: 1,
			'cat'	: 'II/293/glimpse',
			'rad'	: 5,
			'mrk'	: {'m': '>', 'c':'k', 's':9},
			'flt'	: ['3.6mag','4.5mag','5.8mag','8.0mag','e_3.6mag','e_4.5mag','e_5.8mag','e_8.0mag']
			},
	'AKARI':	{# AKARI/IRC mid-IR all-sky Survey (ISAS/JAXA, 2010) 
			'act'	: 1,		
			'fit'	: 1,
			'cat'	: 'II/297/irc',
			'rad'	: 5,
			'mrk'	: {'m': 's', 'c':'m', 's':9},
			'flt'	: ['S09','S18','e_S09','e_S18']
			},
	'MSX6C':	{
			'act'	: 1,		
			'fit'	: 1,
			'cat'	: 'V/114/msx6_gp',
			'rad'	: 5,
			'mrk'	: {'m': 'D', 'c':'b', 's':8},
			'flt'	: ['B1','B2','A','C','D','E','e_B1','e_B2','e_A','e_C','e_D','e_E']
			},
	'WISE':		{# AllWISE Data Release (Cutri+ 2013)
			'act'	: 1,		
			'fit'	: 1,
			'cat'	: 'II/328/allwise',
			'rad'	: 4,
			'mrk'	: {'m': 'd', 'c':'g', 's':9},
			'flt'	: ['W1mag','W2mag','W3mag','W4mag','e_W1mag','e_W2mag','e_W3mag','e_W4mag']
			},
	'IRAS':		{# IRAS catalogue of Point Sources, Version 2.0 (IPAC 1986) 
			'act'	: 1,		
			'fit'	: 0,
			'cat'	: 'II/125/main',
			'rad'	: 4,
			'mrk'	: {'m': '*', 'c':'r', 's':9},
			'flt'	: ['Fnu_12','e_Fnu_12','Fnu_25','e_Fnu_25','Fnu_60','e_Fnu_60','Fnu_100','e_Fnu_100']
			},
		}

#Effective wavelengths (A) of filters and zero point (ZPv) (Source: VOSA)
FLZ_DIR =	{
			'Umag': {'l': 3551.05, 'z': 1438.73},
			'Bmag': {'l': 4369.53, 'z': 4323.91},
			'Vmag': {'l': 5467.57, 'z': 3617.50},
			'Rmag': {'l': 6695.83, 'z': 2908.68},
			'gmag': {'l': 4810.16, 'z': 3631},	#Pan-STARRS is in AB system
			'rmag': {'l': 6155.47, 'z': 3631},	#Pan-STARRS is in AB system
			'imag': {'l': 7503.03, 'z': 3631},	#Pan-STARRS is in AB system
			'zmag': {'l': 8668.36, 'z': 3631},	#Pan-STARRS is in AB system
			'ymag': {'l': 9613.60, 'z': 3631},	#Pan-STARRS is in AB system
			'gpmag': {'l': 4723.59, 'z': 3631},	#Sloan is in AB system
			'rpmag': {'l': 6201.71, 'z': 3631},	#Sloan is in AB system
			'ipmag': {'l': 7672.59, 'z': 3631},	#Sloan is in AB system
			'BPmag': {'l': 5035.75, 'z': 3552.01},
			'Gmag': {'l': 5822.39, 'z': 3228.75},
			'RPmag': {'l': 7619.96, 'z': 2554.95},
			'phot_g_mean_mag': {'l': 5050, 'z': 3534.74},
			'phot_bp_mean_mag': {'l': 6230, 'z': 3296.20},
			'phot_rp_mean_mag': {'l': 7730, 'z': 2620.25},
			'Fmag': {'l': 4977.08, 'z': 3893.51},
			'Nmag': {'l': 6698.86, 'z': 2899.02},
			'Jmag': {'l': 12350, 'z': 1594},
			'Hmag': {'l': 16620, 'z': 1024},
			'Kmag': {'l': 21590, 'z': 666.8},
			'3.6mag': {'l': 35074.83, 'z': 274.53},
			'4.5mag': {'l': 44365.56, 'z': 177.66},
			'5.8mag': {'l': 56280.62, 'z': 113.56},
			'8.0mag': {'l': 75890.54, 'z': 63.70},
			'S09': {'l': 82281.47, 'z': 53.49},	
			'S18': {'l': 176079.23, 'z': 12.05},
			'B1': {'l': 42922.03, 'z': 190.80},
			'B2': {'l': 43492.80, 'z': 185.13},
			'A': {'l': 79512.25, 'z': 57.36},
			'C': {'l': 120723.41, 'z': 26.00},
			'D': {'l': 145915.34, 'z': 17.93},
			'E': {'l': 210177.06, 'z': 8.64},	
			'W1mag': {'l': 33526, 'z': 309.54},
			'W2mag': {'l': 46028, 'z': 171.79},
			'W3mag': {'l': 115608, 'z': 31.67},
			'W4mag': {'l': 220883, 'z': 8.36},
			'Fnu_12': {'l': 101639.68, 'z': 35.36},
			'Fnu_25': {'l': 217351.97, 'z': 7.95},
			'Fnu_60': {'l': 521374.18, 'z': 1.34},
			'Fnu_100': {'l': 954836.89, 'z': 0.406},
		}

GR_SED = 	{ 
		1 : {'nrows': 1, 'ncols': 1, 'figsize': (12,8)},
		2 : {'nrows': 1, 'ncols': 2, 'figsize': (14,6)},
		3 : {'nrows': 1, 'ncols': 3, 'figsize': (21,6)},
		4 : {'nrows': 2, 'ncols': 2, 'figsize': (12,12)},
		5 : {'nrows': 2, 'ncols': 3, 'figsize': (12,10)},
		6 : {'nrows': 2, 'ncols': 3, 'figsize': (15,10)},
		7 : {'nrows': 3, 'ncols': 3, 'figsize': (12,12)},
		8 : {'nrows': 3, 'ncols': 3, 'figsize': (12,12)},
		9 : {'nrows': 3, 'ncols': 3, 'figsize': (15,12)},
		}

MODEL_LEG =	{
		'kurucz' : 'ATLAS9',
		'powr'	 : 'PoWR'
		}

SED_BODIES =	{
		'psph'	:	{
				 'rad':{'name':'rad','value':15,'min':0,'max':400,'vary':True},
				 'ext':{'name':'ext','value': 1.,'min':0,'max':7,'vary':True},
				'teff':{'name':'teff','value':15000,'vary':False},
				'dist':{'name':'dist','value':15000,'vary':False}
				},
		'hd'	:	{
				 'hds':{'name':'hds','value':0.,'vary':True,'min':0.},
				 'hdt':{'name':'hdt','value':1000,'min':600,'max':1500,'vary':True}
				},
		'wd'	:	{
				 'wds':{'name':'wds','value':0.,'vary':True,'min':0.},
				 'wdt':{'name':'wdt','value':500,'min':150,'max':600,'vary':True}
				}
		}	
