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
		{'legend.fontsize': 12,
	 	'font.size':  12,
         	'axes.labelsize': 16,
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
		{'legend.fontsize': 10,
	 	'font.size':  8,
         	'axes.labelsize': 14,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 18,
         	'ytick.labelsize': 11},
'cr'		:
		{'legend.fontsize': 11,
		#'text.usetex' : True,
	 	'font.size':  11,
         	'axes.labelsize': 16,
         	'axes.titlesize': 13,
         	'xtick.labelsize': 15,
         	'ytick.labelsize': 15},
'panel'		:
		{'legend.fontsize': 12,
	 	'font.size':  14,
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
		'sed': r'$f_{\lambda}$ (erg/s/cm$^2$/A)'
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
		'LOGG' : r'log$g$', 'MASS' : '$M_{evol}$', 'VMAC': r'$v_{mac}$ [km s$^{-1}$]', 'VMIC' : r'$v_{mic}$ [km s$^{-1}$]',
		'NABUN' : r'N/H', 'LOGD' : r'log$D$', 'S_MASS' : r'$M_{Rg}$',
		'LOGW' : r'log($W$ [mag])', 'LOGR0' : r'log($R_{0}$ [mag])', 'TAU' : r'$\tau$ [d]', 'GAMMA' : r'$\gamma$' ,
		'TEFF' :  r'log(T$_{\rm eff}$ [K])', 'TESS_time' : r'Time $-$ 2457000 [BTJD d]',
		'TESS_freq' : r'Frequency [d$^{-1}$]','S_LOGL' : r'log($L$/L$_{\odot})$',
		'MAD' : r'MAD', 'SVAR': r'$\sigma$ [mag]', 'ZCROSS' : r'$D_{0}$', 'PSI': r'$\psi^2$',
		'ETA' : r'log($\eta$)', 'SKEW': r'skw', 'A_V' : r'$A_{V}$', 'EDD' : r'$\Gamma_{e}$',
		'GABS' : r'$M_{G}$ [mag]', 'VCHAR' : r'log($\nu_{char}$ [d$^{-1}$])',
		'FF' : r'$f_{i}$ [d$^{-1}$]', 'A_FF' : r'$A_{i}$ [mag]', 'HF' : r'$jf$', 'FFR' : r'$f_{1}/f_{2}$', 
		'A_FFR' : r'$A_{f_{1}}/A_{f_{2}}$', 'BETA' : r'$\beta$', 'VINF' : r'$v_{inf}$ [km s$^{-1}$]',
		'INDFF' : r'#$f_{i}$','INDFFS' : r'#$f_{i,sec}$'}

FMT_PR = 	{
		'TEFF' : '%.2f','VSINI' : '%d','MDOT' : '%.2f',
		'LOGLM': '%.2f', 'S_LOGL' : '%.1f', 'LOGQ' : '%.2f', 
		'LOGG' : '%.1f', 'MASS' : '%d', 'VMAC': '%d', 'VMIC' : '%d', 
		'NABUN' : '%.1f', 'LOGD' : '%.2f', 'LUM': '%.1f', 'S_MASS' : '%d',
		'EDD' : '%.2f',
		}


BR_AX_C = .02

LSPRW_PLOT_NUM = 4

HIST_COLORS = ['dodgerblue','limegreen','r','c','k','b','m']

BIN_PROP =	{
		'LOW' : 0.05, 
		'UP' : 1.31, 
		'RES' : 0.05
		}

CBAR_TICK_SIZE = 17
CBAR_TITLE_SIZE = 18

HR_MS = 14
HR_YLIM = [4.21, 6.45]
HR_XLIM = [4.78,3.51]

SHR_MS = 14
SHR_YLIM = [4.55,0.8]
SHR_XLIM = [4.77,3.97]

SED_YLIM = [8e-18,9e-8]
SED_XLIM = [7e+2,4e+5]

HD_LMT = 	{ 
		'TEFF' : [3.4,3.8,4.62],
	     	'LOGL' : [5.8,5.8,6.5],
		'lw'   : 3,
		'c'    : 'k'	
	   	}

PRSG_AGE =	{
		'20'	: 9.86029e+06,
		'25'	: 7.99503e+06,
		'32'	: 6.70746e+06,
		'40'	: 5.74518e+06,
		'60'	: 4.49777e+06,
		'120'	: 3.05433e+06
		}

