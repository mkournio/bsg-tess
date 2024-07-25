from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from constants.c_plot import *
import os


class GridTemplate(object):

	# Class for the creation and management of plotting grid

	def __init__(self, rows_page = 3, cols_page = 1, output_format = 'pdf', params = PLOT_PARAMS['lc'], inter = False, figsize = SIZE_GRID, **kwargs):

		plt.rcParams.update(params)

        	self.__dict__.update(kwargs)
		self.rows_page = rows_page
		self.cols_page = cols_page
		self.output_format = output_format
		self.inter = inter
		self.figsize = figsize

		self.filename = kwargs.pop('output_name','GridPlot')
		self.fig_xlabel  = kwargs.pop('fig_xlabel','X LABEL')
		self.fig_ylabel  = kwargs.pop('fig_ylabel','Y LABEL')
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
		self.ax_pos = plot_row

		if plot_row < self.rows_page - 1:
			if self.coll_x:
			 	ax.tick_params(labelbottom=False)
			 	self.gs.update(hspace=0.05)
		else:
			if 'col_labels' in self.__dict__:
				ax.set_xlabel(STY_LB[self.col_labels[plot_col]])

		if plot_row == 0 and 'sup_xlabels' in self.__dict__:
				ax.set_title(self.sup_xlabels[plot_col])

		if plot_col > 0 :
			if self.coll_y:
				ax.tick_params(labelleft=False)
			 	self.gs.update(wspace=0.05)
		else:
			if 'row_labels' in self.__dict__:
				ax.set_ylabel(STY_LB[self.row_labels[plot_row]])			

		if 'xlim' in self.__dict__:ax.set_xlim(self.xlim)
		if 'ylim' in self.__dict__:ax.set_ylim(self.ylim)

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
		while os.path.exists('%s_%s.%s' % (self.filename,fl_id,self.output_format)): fl_id += 1

		return '%s_%s' % (self.filename,fl_id)		

	def _grid_new(self):
		
		#Create grid on new page

		self.fig = plt.figure(figsize=self.figsize)
		self.gs = GridSpec(self.rows_page, self.cols_page, figure=self.fig)

		return

	def _page_close(self):

		#Close/Save grid in page
		self._save_output()
		if self.inter: plt.show()
		self.fig.clf(); plt.close(self.fig)

		return

	def _save_output(self):

		self.glob_ax = self.fig.add_subplot(self.gs[:self.ax_pos+1, :], frameon=False)
		self.glob_ax.tick_params(labelleft = False, labelbottom = False, bottom = False, left = False)
		self.glob_ax.set_xlabel(self.fig_xlabel, fontsize = SIZE_XLABEL_FIG, labelpad=30)
		self.glob_ax.set_ylabel(self.fig_ylabel, fontsize = SIZE_YLABEL_FIG, labelpad=55)

		if self.output_format != None:
			filename = '%s.%s' % (self._get_filename(),self.output_format)
			self.fig.savefig('%s.%s' % (self._get_filename(),self.output_format), format = self.output_format,bbox_inches='tight')

		return					


