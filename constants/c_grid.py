#### Plotting constants and labels
grid_params = {'legend.fontsize': 8,
	 	'font.size':  15,
         	'axes.labelsize': 12,
         	'axes.titlesize': 11,
         	'xtick.labelsize': 12,
         	'ytick.labelsize': 9}

PLOT_XLC_NCOL = 3

PLOT_XLABEL =   { 
		'lc' : r'Time $-$ 2457000 [BTJD d]',
		'ls' : r'Frequency [d$^{-1}$]'
		}
PLOT_YLABEL =   { 
		'lc' : r'mag',
		'ls' : r'Power'
		}
LC_COLOR = 	{
		'SPOC'	 : 'c',
		'FFI'	 : 'k',
		'polyfit': 'r'
		}
AX_FMARK =	{
		'(F)' : {'xmin' : 0.7, 'xmax' : 0.9, 'ls' : '-', 'lw' : 3},
		'(H)' : {'xmin' : 0.7, 'xmax' : 0.9, 'ls' : '--', 'lw' : 1.5},
		'(C)' : {'xmin' : 0.7, 'xmax' : 0.9, 'ls' : '-', 'lw' : 0.2}
		}
BR_AX_C = .02
LSPRW_PLOT_NUM = 4


