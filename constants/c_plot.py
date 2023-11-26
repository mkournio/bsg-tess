#### Plotting constants and labels
grid_params = {'legend.fontsize': 8,
	 	'font.size':  11,
         	'axes.labelsize': 14,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 12,
         	'ytick.labelsize': 9}

panel_params = {'legend.fontsize': 8,
	 	'font.size':  11,
         	'axes.labelsize': 16,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 16,
         	'ytick.labelsize': 16}


PLOT_XLC_NCOL = 3

PLOT_XLABEL =   { 
		'lc' : r'Time $-$ 2457000 [BTJD d]',
		'ls' : r'Frequency [d$^{-1}$]',
		'sed': r'Wavelength (A)'
		}
PLOT_YLABEL =   { 
		'lc' : r'mag',
		'ls' : r'Power',
		'sed': r'$F_{\lambda}$ (erg/s/cm$^2$/A)'
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

STY_LB = 	{
		'TEFF' : r'T$_{\rm eff}$','VSINI' : r'$v$sin $i$','MDOT' : r'$\dot{M}$',
		'LOGLM': r'log(L/M)', 'LOGL' : r'log(L/L$_{\odot}$)', 'LOGQ' : r'log$Q$',
		'LOGG' : r'log$g$', 'MASS' : 'M', 'VMAC': r'$v_{mac}$', 'VMIC' : r'$v_{mic}$',
		'NABUN' : r'N/H', 'LOGD' : r'log$D$', 'LUM': r'log(L/L$_{\odot})_{\rm SED}$', 
		'w' : r'$W$', 'zero' : r'$R_{0}$', 'tau' : r'$\tau$', 'gamma' : r'$\gamma$' ,
		'LTEFF' :  r'log(T$_{\rm eff}$ [K])'
		}

FMT_PR = 	{
		'TEFF' : '%d','VSINI' : '%d','MDOT' : '%.2f',
		'LOGLM': '%.2f', 'LOGL' : '%.1f', 'LOGQ' : '%.2f', 
		'LOGG' : '%.1f', 'MASS' : '%d', 'VMAC': '%d', 'VMIC' : '%d', 
		'NABUN' : '%.1f', 'LOGD' : '%.2f', 'LUM': '%.1f' 
		}


BR_AX_C = .02

LSPRW_PLOT_NUM = 4

HIST_COLORS = ['r','limegreen','dodgerblue']

BIN_PROP =	{
		'LOW' : 0.05, 
		'UP' : 1.50, 
		'RES' : 0.05
		}

HD_LMT = 	{ 
		'TEFF' : [3.5,3.8,4.45],
	     	'LOGL' : [5.8,5.8,6.35],
		'lw'   : 3,
		'c'    : 'k'	
	   	}

