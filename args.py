filename = 'LaPlata_input'

xmtabs = ['HAUCKE+19','FRASER+10']
xmcols = ['TEFF','NABUN','MASS','LOGL','VMIC', 'LOGG','VMAC','MDOT']

# CURRENTLY VIZIER TABLES HAVE TO ME IMPORTED ONE BY ONE
xmviz = {
	'vizier:I/352/gedr3dis' : {'xcols' : 'DIST', 'viz_col' : 'rpgeo'},
	'vizier:IV/38/tic'      : {'xcols' : 'TIC', 'viz_col' : 'TIC'},
	'vizier:II/328/allwise' : {'xcols' : 'WISEFL', 'viz_col' : 'qph', 'column_format' : 'U5'},
	'vizier:IV/39/tic82'    : {'xcols' : 'TMAG', 'viz_col' : 'Tmag'},
	'vizier:I/311/hip2'     : {'xcols' : 'HDIST', 'viz_col' : 'Plx'}
	   }
	
