TEX_SAMPLE_TAB=	{ 
		'preamble': r'\centering',
	        'tabletype': 'table*',
		'data_end': r'\hline',
		'header_end': r'\hline', 
		'header_start': r'\hline\hline',
		}

TEX_FREQ_TAB = {
		'PRE'  :  ("\\begin{table*}\n"
			   "\\caption{\\label{sample} List of frequencies.}\n"
			   "\\begin{tabular}{l c c c c}\n"
			   "\\hline\\\n"
			   "Star & Sector & ID & Frequency & SNR\\\\\n"
			   "\\hline\n") ,
		'POST' :  ("\\hline\n"
			   "\\end{tabular}\n"
			   "\\end{table*}")
		}
