from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from constants.c_grid import *
import os

plt.rcParams.update(grid_params)

class GridTemplate(object):

	# Class for the creation and management of plotting grid

	def __init__(self, rows_page = 3, cols_page = 1, output_format = 'pdf', inter = False, **kwargs):

        	self.__dict__.update(kwargs)
		self.rows_page = rows_page
		self.cols_page = cols_page
		self.output_format = output_format
		self.inter = inter

		self.filename = kwargs.pop('output_name','GridPlot')
		self.coll_x = kwargs.pop('coll_x',False)
		self.coll_y = kwargs.pop('coll_y',False)

		#Initiating pdf
		if self.output_format == 'pdf': 
				self.pdf = PdfPages(self._get_filename())

		self.ind = 0

	def GridAx(self):

		# Returns the position axis on the grid when called from the child class

		plot_row = (self.ind / self.cols_page) % self.rows_page
		plot_col = self.ind % self.cols_page
		plot_tot = (self.cols_page*self.rows_page)

		if (self.ind % plot_tot  == 0) :
			if self.ind > 0 : self._page_close()	# Close/Save the existing grid
			self._grid_new()

		ax = self.fig.add_subplot(self.gs[plot_row,plot_col])


		if plot_row < self.rows_page - 1:
			if self.coll_x:
			 	ax.tick_params(labelbottom=False)
			 	self.gs.update(hspace=0.1)
		else:
			if 'col_labels' in self.__dict__:
				ax.set_xlabel(STY_LB[self.col_labels[plot_col]])

		if plot_row == 0 and 'sup_xlabels' in self.__dict__:
				ax.set_title(self.sup_xlabels[plot_col])

		if plot_col > 0 :
			if self.coll_y:
				ax.tick_params(labelleft=False)
			 	self.gs.update(wspace=0.1)
		else:
			if 'row_labels' in self.__dict__:
				ax.set_ylabel(STY_LB[self.row_labels[plot_row]])
			

		self.ind += 1

		return ax



	def GridClose(self):
		
		self._save_output()
		if self.output_format == 'pdf': self.pdf.close()
		if self.inter: plt.show()
		self.fig.clf(); plt.close('all')

		return

	def _get_filename(self):

		fl_id = 0 
		while os.path.exists('%s_%s' % (self.filename,fl_id)): fl_id += 1

		return '%s_%s' % (self.filename,fl_id)		

	def _grid_new(self):
		
		#Create grid on new page

		self.fig = plt.figure(figsize=(16,10))
		self.gs = GridSpec(self.rows_page, self.cols_page, figure=self.fig)
		if 'fig_xlabel' in self.__dict__:
			self.fig.text(0.5,0.04,self.fig_xlabel,size = 15,ha="center",va="center")
		if 'fig_ylabel' in self.__dict__:
			self.fig.text(0.06,0.5,self.fig_ylabel,size = 15, ha="center", va="center", rotation=90)

		return

	def _page_close(self):

		#Close/Save grid in page
		self._save_output()
		if self.inter: plt.show()
		self.fig.clf(); plt.close(self.fig)

		return

	def _save_output(self):

		if self.output_format == 'pdf':
			self.pdf.savefig(self.fig)
		elif self.output_format == 'eps':
			self.fig.savefig('%s' % self._get_filename(), format = 'eps')	

		return					


