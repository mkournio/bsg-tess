

def colorbar(fig,vmin,vmax,label,cbar='jet'):

		import matplotlib.cm as cm
		from matplotlib.colors import Normalize

		cmap = cm.get_cmap(cbar).reversed()
		cmap.set_under('white')
		norm = Normalize(vmin,vmax)

		cbar_ax = fig.add_axes([0.91, 0.11, 0.03, 0.77])
		cmmapable = cm.ScalarMappable(norm, cmap)
		cmmapable.set_array([])
		cbar=fig.colorbar(cmmapable, cax=cbar_ax)
		cbar.ax.set_title(label)
		cbar.ax.tick_params(labelsize=14)

		return cmap, norm
