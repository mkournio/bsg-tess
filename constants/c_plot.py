#### Plotting constants and labels
SIZE_FONT_SUB = 12
SIZE_XLABEL_FIG = 22
SIZE_YLABEL_FIG = 22
SIZE_XLABEL_SUB = 11
SIZE_YLABEL_SUB = 11

SIZE_GRID = (16,20)

PLOT_PARAMS =	{
'lc'		:
		{'legend.fontsize': 8,
	 	'font.size':  SIZE_FONT_SUB,
         	'axes.labelsize': 14,
         	'axes.titlesize': 13,
         	'xtick.labelsize': SIZE_XLABEL_SUB,
         	'ytick.labelsize': SIZE_YLABEL_SUB},
'ls'		:
		{'legend.fontsize': 8,
	 	'font.size':  13,
         	'axes.labelsize': 14,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 17,
         	'ytick.labelsize': 17},
'prew'		:
		{'legend.fontsize': 8,
	 	'font.size':  SIZE_FONT_SUB+6,
         	'axes.labelsize': 24,
         	'axes.titlesize': 13,
         	'xtick.labelsize': SIZE_XLABEL_SUB+5,
         	'ytick.labelsize': SIZE_YLABEL_SUB+5},
'sed'		:
		{'legend.fontsize': 8,
	 	'font.size':  10,
         	'axes.labelsize': 14,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 18,
         	'ytick.labelsize': 11},
'cr'		:
		{'legend.fontsize': 11,
	 	'font.size':  13,
         	'axes.labelsize': 16,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 15,
         	'ytick.labelsize': 15},
'panel'		:
		{'legend.fontsize': 8,
	 	'font.size':  11,
         	'axes.labelsize': 22,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 23,
         	'ytick.labelsize': 23}
		}

PLOT_XLC_NCOL = 3

PLOT_XLABEL =   { 
		'lc' : r'Time $-$ 2457000 [BTJD d]',
		'ls' : r'Frequency [d$^{-1}$]',
		'sed': r'Wavelength (A)'
		}

PLOT_YLABEL =   { 
		'lc' : r'$\Delta$m (mag)',
		'ls' : r'Amplitude (mag)',
		'sed': r'$F_{\lambda}$ (erg/s/cm$^2$/A)'
		}



LC_COLOR = 	{
		'SPOC'	 : 'c',
		'FFI'	 : 'r',
		'any': 'k'
		}
AX_FMARK =	{
		'(F)' : {'xmin' : 0.7, 'xmax' : 0.9, 'ls' : '-', 'lw' : 3},
		'(H)' : {'xmin' : 0.7, 'xmax' : 0.9, 'ls' : '--', 'lw' : 1.5},
		'(C)' : {'xmin' : 0.7, 'xmax' : 0.9, 'ls' : '-', 'lw' : 0.2}
		}

STY_LB = 	{
		'VSINI' : r'$v$sin $i$','MDOT' : r'$\dot{M}$',
		'LOGLM': r'log(L/M)', 'LOGL' : r'log(L/L$_{\odot}$)', 'LOGQ' : r'log$Q$',
		'LOGG' : r'log$g$', 'MASS' : '$M_{evol}$', 'VMAC': r'$v_{mac}$', 'VMIC' : r'$v_{mic}$',
		'NABUN' : r'N/H', 'LOGD' : r'log$D$', 'S_MASS' : r'$M_{Rg}$',
		'w' : r'log($W$ [mag])', 'zero' : r'log($R_{0}$ [mag])', 'tau' : r'log($\tau$ [d])', 'gamma' : r'$\gamma$' ,
		'TEFF' :  r'log(T$_{\rm eff}$ [K])', 'TESS_time' : r'Time $-$ 2457000 [BTJD d]',
		'TESS_freq' : r'Frequency [d$^{-1}$]','S_LOGL' : r'log($L$/L$_{\odot})_{Gaia}$',
		'MAD' : r'MAD', 'SVAR': r'$\sigma$', 'ZCROSS' : r'$D_{0}$', 'PSI': r'$\psi^2$',
		'ETA' : r'log($\eta$)', 'SKEW': r'skw', 'A_V' : r'$A_{V}$', 'EDD' : r'$\Gamma_{e}$',
		'GABS' : r'$M_{G}$ [mag]', 'VCHAR' : r'log($\nu_{char}$ [d$^{-1}$])',
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

CBAR_TICK_SIZE = 17
CBAR_TITLE_SIZE = 18

HR_YLIM = [4.21, 6.45]

HR_XLIM = [4.78,3.51]

HD_LMT = 	{ 
		'TEFF' : [3.4,3.8,4.62],
	     	'LOGL' : [5.8,5.8,6.5],
		'lw'   : 3,
		'c'    : 'k'	
	   	}

