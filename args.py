filename = 'LaPlata'

xmtabs = ['MATTIAS','HAUCKE+19','FRASER+10']
xmcols = ['TEFF','NABUN','MASS','LOGL','VMIC', 'LOGG','VMAC','VSINI','MDOT','BETA','VINF']

# CURRENTLY VIZIER TABLES HAVE TO ME IMPORTED ONE BY ONE
xmviz = {
	'GAIA_1'		: {'cat' : 'vizier:I/352/gedr3dis', 'xcols' : 'DIST', 'viz_col' : 'rpgeo'},
	'GAIA_2'		: {'cat' : 'vizier:I/352/gedr3dis', 'xcols' : 'bDIST', 'viz_col' : 'b_rpgeo'},
	'GAIA_3'		: {'cat' : 'vizier:I/352/gedr3dis', 'xcols' : 'BDIST', 'viz_col' : 'B_rpgeo'},
	'TESS'			: {'cat' : 'vizier:IV/38/tic', 'xcols' : 'TIC', 'viz_col' : 'TIC'},
	'WISE'			: {'cat' : 'vizier:II/328/allwise', 'xcols' : 'WISEFL', 'viz_col' : 'qph', 'column_format' : 'U5'},
	'TESS_2'		: {'cat' : 'vizier:IV/39/tic82', 'xcols' : 'TMAG', 'viz_col' : 'Tmag'},
	'HIPPARCOS'		: {'cat' : 'vizier:I/311/hip2', 'xcols' : 'HDIST', 'viz_col' : 'Plx'},
	'HIPPARCOS_2'		: {'cat' : 'vizier:I/311/hip2', 'xcols' : 'e_HDIST', 'viz_col' : 'e_Plx'}
	   }
	
