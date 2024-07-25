from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from constants.c_plot import *
import os



class PanelTemplate(object):

	# Class for the creation and management of plotting grid
	def __init__(self, output_format = None, params = PLOT_PARAMS['panel'],  inter = False, **kwargs):

		plt.rcParams.update(params)

        	self.__dict__.update(kwargs)
		self.output_format = output_format
		self.inter = inter

		self.inv_x = kwargs.pop('inv_x',False)
		self.inv_y = kwargs.pop('inv_y',False)

		self.filename = kwargs.pop('output_name','PanelPlot')

		#Initiating pdf
		if self.output_format == 'pdf': 
				self.pdf = PdfPages(self._get_filename())

	def PanelAx(self):

		self.fig, ax = plt.subplots(figsize=(10, 10))

		if 'x_label' in self.__dict__: ax.set_xlabel(STY_LB[self.x_label])
		if 'y_label' in self.__dict__: ax.set_ylabel(STY_LB[self.y_label])
		
		if self.inv_x : ax.invert_xaxis()
		if self.inv_y : ax.invert_yaxis()

		return ax

	def PanelClose(self):
		
		self._save_output()
		if self.output_format == 'pdf': self.pdf.close()
		if self.inter: plt.show()
		self.fig.clf(); plt.close('all')

		return

	def _get_filename(self):

		fl_id = 0 
		while os.path.exists('%s_%s' % (self.filename,fl_id)): fl_id += 1

		return '%s_%s' % (self.filename,fl_id)		

	def _save_output(self):

		filename = '%s.%s' % (self._get_filename(),self.output_format)
		self.fig.savefig(filename, format = self.output_format, bbox_inches='tight')	

		return					





