from __future__ import division, print_function
from abc import ABCMeta, abstractmethod
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib.ticker import MaxNLocator, Formatter, Locator
from matplotlib.widgets import Slider, Button
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from tadtool.tad import GenomicRegion, sub_matrix_regions, sub_data_regions, \
    data_array, insulation_index, sub_vector_regions, sub_regions, \
    call_tads_insulation_index, directionality_index, call_tads_directionality_index, normalised_insulation_index
import math
import copy
import numpy as np
from bisect import bisect_left
from future.utils import string_types

try:
    import Tkinter as tk
    import tkFileDialog as filedialog
except ImportError:
    import tkinter as tk
    from tkinter import filedialog


class BasePlotter(object):

    __metaclass__ = ABCMeta

    def __init__(self, title):
        self._ax = None
        self.cax = None
        self.title = title

    @abstractmethod
    def _plot(self, region=None, **kwargs):
        raise NotImplementedError("Subclasses need to override _plot function")

    @abstractmethod
    def plot(self, region=None, **kwargs):
        raise NotImplementedError("Subclasses need to override plot function")

    @property
    def fig(self):
        return self._ax.figure

    @property
    def ax(self):
        if not self._ax:
            _, self._ax = plt.subplots()
        return self._ax

    @ax.setter
    def ax(self, value):
        self._ax = value


class GenomeCoordFormatter(Formatter):
    """
    Process axis tick labels to give nice representations
    of genomic coordinates
    """
    def __init__(self, chromosome, display_scale=True):
        """
        :param chromosome: :class:`~kaic.data.genomic.GenomicRegion` or string
        :param display_scale: Boolean
                              Display distance scale at bottom right
        """
        if isinstance(chromosome, GenomicRegion):
            self.chromosome = chromosome.chromosome
        else:
            self.chromosome = chromosome
        self.display_scale = display_scale

    def _format_val(self, x, prec_offset=0):
        if x == 0:
            oom_loc = 0
        else:
            oom_loc = int(math.floor(math.log10(abs(x))))
        view_range = self.axis.axes.get_xlim()
        oom_range = int(math.floor(math.log10(abs(view_range[1] - view_range[0]))))
        if oom_loc >= 3:
            return "{:.{prec}f}kb".format(x/1000, prec=max(0, 3 + prec_offset - oom_range))
        return "{:.0f}b".format(x)

    def __call__(self, x, pos=None):
        """
        Return label for tick at coordinate x. Relative position of
        ticks can be specified with pos. First tick gets chromosome name.
        """
        s = self._format_val(x, prec_offset=1)
        if pos == 0 or x == 0:
            return "{}:{}".format(self.chromosome, s)
        return s

    def get_offset(self):
        """
        Return information about the distances between
        tick bars and the size of the view window.
        Is called by matplotlib and displayed in lower right corner
        of plots.
        """
        if not self.display_scale:
            return ""
        view_range = self.axis.axes.get_xlim()
        view_dist = abs(view_range[1] - view_range[0])
        tick_dist = self.locs[2] - self.locs[1]
        minor_tick_dist = tick_dist/5
        minor_tick_dist_str = self._format_val(minor_tick_dist, prec_offset=2)
        tick_dist_str = self._format_val(tick_dist, prec_offset=1)
        view_dist_str = self._format_val(view_dist)
        return "{}|{}|{}".format(minor_tick_dist_str, tick_dist_str, view_dist_str)


class GenomeCoordLocator(MaxNLocator):
    """
    Choose locations of genomic coordinate ticks on the plot axis.
    Behaves like default Matplotlib locator, except that it always
    places a tick at the start and the end of the window.
    """
    def __call__(self):
        vmin, vmax = self.axis.get_view_interval()
        ticks = self.tick_values(vmin, vmax)
        # Make sure that first and last tick are the start
        # and the end of the genomic range plotted. If next
        # ticks are too close, remove them.
        ticks[0] = vmin
        ticks[-1] = vmax
        if ticks[1] - vmin < (vmax - vmin)/(self._nbins*3):
            ticks = np.delete(ticks, 1)
        if vmax - ticks[-2] < (vmax - vmin)/(self._nbins*3):
            ticks = np.delete(ticks, -2)
        return self.raise_if_exceeds(np.array(ticks))


class MinorGenomeCoordLocator(Locator):
    """
    Choose locations of minor tick marks between major
    tick labels. Modification of the Matplotlib AutoMinorLocator,
    except that it uses the distance between 2nd and 3rd major
    mark as reference, instead of 2nd and 3rd.
    """
    def __init__(self, n):
        self.ndivs = n

    def __call__(self):
        majorlocs = self.axis.get_majorticklocs()
        try:
            majorstep = majorlocs[2] - majorlocs[1]
        except IndexError:
            # Need at least two major ticks to find minor tick locations
            # TODO: Figure out a way to still be able to display minor
            # ticks without two major ticks visible. For now, just display
            # no ticks at all.
            majorstep = 0
        if self.ndivs is None:
            if majorstep == 0:
                # TODO: Need a better way to figure out ndivs
                ndivs = 1
            else:
                x = int(np.round(10 ** (np.log10(majorstep) % 1)))
                if x in [1, 5, 10]:
                    ndivs = 5
                else:
                    ndivs = 4
        else:
            ndivs = self.ndivs
        minorstep = majorstep / ndivs
        vmin, vmax = self.axis.get_view_interval()
        if vmin > vmax:
            vmin, vmax = vmax, vmin
        if len(majorlocs) > 0:
            t0 = majorlocs[1]
            tmin = ((vmin - t0) // minorstep + 1) * minorstep
            tmax = ((vmax - t0) // minorstep + 1) * minorstep
            locs = np.arange(tmin, tmax, minorstep) + t0
            cond = np.abs((locs - t0) % majorstep) > minorstep / 10.0
            locs = locs.compress(cond)
        else:
            locs = []
        return self.raise_if_exceeds(np.array(locs))


class BasePlotter1D(BasePlotter):

    __metaclass__ = ABCMeta

    def __init__(self, title):
        BasePlotter.__init__(self, title=title)

    def plot(self, region=None, ax=None, **kwargs):
        if isinstance(region, string_types):
            region = GenomicRegion.from_string(region)
        if ax:
            self.ax = ax
        # set genome tick formatter
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(region))
        self.ax.xaxis.set_major_locator(GenomeCoordLocator(nbins=5))
        self.ax.xaxis.set_minor_locator(MinorGenomeCoordLocator(n=5))
        self.ax.set_title(self.title)
        self._plot(region, **kwargs)
        self.ax.set_xlim(region.start, region.end)
        return self.fig, self.ax


def prepare_normalization(norm="lin", vmin=None, vmax=None):
    if isinstance(norm, mpl.colors.Normalize):
        norm.vmin = vmin
        norm.vmax = vmax
        return norm
    if norm == "log":
        return mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    elif norm == "lin":
        return mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("'{}'' not a valid normalization method.".format(norm))


class BasePlotterHic(object):

    __metaclass__ = ABCMeta

    def __init__(self, hic_matrix, regions=None, colormap='RdBu', norm="log",
                 vmin=None, vmax=None, show_colorbar=True, blend_masked=False):
        if regions is None:
            for i in range(hic_matrix.shape[0]):
                regions.append(GenomicRegion(chromosome='', start=i, end=i))
        self.regions = regions
        self.hic_matrix = hic_matrix

        self.colormap = copy.copy(mpl.cm.get_cmap(colormap))
        if blend_masked:
            self.colormap.set_bad(self.colormap(0))
        self._vmin = vmin
        self._vmax = vmax
        self.norm = prepare_normalization(norm=norm, vmin=vmin, vmax=vmax)
        self.colorbar = None
        self.slider = None
        self.show_colorbar = show_colorbar

    def add_colorbar(self, ax=None):
        ax = self.cax if ax is None else ax
        cmap_data = mpl.cm.ScalarMappable(norm=self.norm, cmap=self.colormap)
        cmap_data.set_array([self.vmin, self.vmax])

        self.colorbar = plt.colorbar(cmap_data, cax=ax, orientation="vertical")

    @property
    def vmin(self):
        return self._vmin if self._vmin else np.nanmin(self.hic_matrix)

    @property
    def vmax(self):
        return self._vmax if self._vmax else np.nanmax(self.hic_matrix)


class HicPlot(BasePlotter1D, BasePlotterHic):
    def __init__(self, hic_matrix, regions=None, title='', colormap='viridis', max_dist=None, norm="log",
                 vmin=None, vmax=None, show_colorbar=True, blend_masked=False):
        BasePlotter1D.__init__(self, title=title)
        BasePlotterHic.__init__(self, hic_matrix, regions=regions, colormap=colormap, vmin=vmin, vmax=vmax,
                                show_colorbar=show_colorbar, norm=norm, blend_masked=blend_masked)
        self.max_dist = max_dist
        self.hicmesh = None

    def _plot(self, region=None, cax=None):
        if region is None:
            raise ValueError("Cannot plot triangle plot for whole genome.")

        hm, sr = sub_matrix_regions(self.hic_matrix, self.regions, region)
        hm[np.tril_indices(hm.shape[0])] = np.nan
        # Remove part of matrix further away than max_dist
        if self.max_dist is not None:
            for i in range(hm.shape[0]):
                i_region = sr[i]
                for j in range(hm.shape[1]):
                    j_region = sr[j]

                    if j_region.start-i_region.end > self.max_dist:
                        hm[i, j] = np.nan

        hm_masked = np.ma.MaskedArray(hm, mask=np.isnan(hm))
        # prepare an array of the corner coordinates of the Hic-matrix
        # Distances have to be scaled by sqrt(2), because the diagonals of the bins
        # are sqrt(2)*len(bin_size)
        sqrt2 = math.sqrt(2)
        bin_coords = np.r_[[(x.start - 1) for x in sr], sr[-1].end]/sqrt2
        X, Y = np.meshgrid(bin_coords, bin_coords)
        # rotatate coordinate matrix 45 degrees
        sin45 = math.sin(math.radians(45))
        X_, Y_ = X*sin45 + Y*sin45, X*sin45 - Y*sin45
        # shift x coords to correct start coordinate and center the diagonal directly on the
        # x-axis
        X_ -= X_[1, 0] - (sr[0].start - 1)
        Y_ -= .5*np.min(Y_) + .5*np.max(Y_)
        # create plot
        self.hicmesh = self.ax.pcolormesh(X_, Y_, hm_masked, cmap=self.colormap, norm=self.norm)

        # set limits and aspect ratio
        #self.ax.set_aspect(aspect="equal")
        ylim_max = 0.5*(region.end-region.start)
        if self.max_dist is not None and self.max_dist/2 < ylim_max:
            ylim_max = self.max_dist/2
        self.ax.set_ylim(0, ylim_max)
        # remove y ticks
        self.ax.set_yticks([])
        # hide background patch
        self.ax.patch.set_visible(False)
        if self.show_colorbar:
            self.add_colorbar(cax)

    def set_clim(self, vmin, vmax):
        self.hicmesh.set_clim(vmin=vmin, vmax=vmax)
        if self.colorbar is not None:
            self.colorbar.vmin = vmin
            self.colorbar.vmax = vmax
            self.colorbar.draw_all()


class DataArrayPlot(BasePlotter1D):
    def __init__(self, data, window_sizes=None, regions=None, title='', midpoint=None,
                 colormap='coolwarm_r', vmax=None, current_window_size=0, log_y=True):
        if regions is None:
            regions = []
            for i in range(data.shape[1]):
                regions.append(GenomicRegion(chromosome='', start=i, end=i))
        self.regions = regions

        BasePlotter1D.__init__(self, title=title)
        self.da = data
        if window_sizes is None:
            window_sizes = []
            try:
                l = len(data)
            except TypeError:
                l = data.shape[0]
            for i in range(l):
                window_sizes.append(i)

        self.window_sizes = window_sizes
        self.colormap = colormap
        self.midpoint = midpoint
        self.mesh = None
        self.vmax = vmax
        self.window_size_line = None
        self.current_window_size = current_window_size
        self.log_y = log_y

    def _plot(self, region=None, cax=None):
        da_sub, regions_sub = sub_data_regions(self.da, self.regions, region)

        da_sub_masked = np.ma.MaskedArray(da_sub, mask=np.isnan(da_sub))
        bin_coords = np.r_[[(x.start - 1) for x in regions_sub[:-1]], regions_sub[-1].end]
        x, y = np.meshgrid(bin_coords, self.window_sizes)
        
        self.mesh = self.ax.pcolormesh(x, y, da_sub_masked, cmap=self.colormap, vmax=self.vmax)
        self.colorbar = plt.colorbar(self.mesh, cax=cax, orientation="vertical")
        self.window_size_line = self.ax.axhline(self.current_window_size, color='red')

        if self.log_y:
            self.ax.set_yscale("log")
        self.ax.set_ylim((np.nanmin(self.window_sizes), np.nanmax(self.window_sizes)))

    def set_clim(self, vmin, vmax):
        self.mesh.set_clim(vmin=vmin, vmax=vmax)
        if self.colorbar is not None:
            self.colorbar.vmin = vmin
            self.colorbar.vmax = vmax
            self.colorbar.draw_all()

    def update(self, window_size):
        self.window_size_line.set_ydata(window_size)


class TADPlot(BasePlotter1D):
    def __init__(self, regions, title='', color='black'):
        BasePlotter1D.__init__(self, title=title)

        self.regions = regions
        self.color = color
        self.current_region = None

    def _plot(self, region=None, cax=None):
        self.current_region = region
        try:
            sr, start_ix, end_ix = sub_regions(self.regions, region)
            trans = self.ax.get_xaxis_transform()
            for r in sr:
                region_patch = patches.Rectangle(
                    (r.start, .2),
                    width=abs(r.end - r.start), height=.6,
                    transform=trans,
                    facecolor=self.color,
                    edgecolor='white',
                    linewidth=2.
                )
                self.ax.add_patch(region_patch)
        except ValueError:
            pass
        self.ax.axis('off')

    def update(self, regions):
        self.regions = regions
        self.ax.cla()
        self.plot(region=self.current_region, ax=self.ax)


class DataLinePlot(BasePlotter1D):
    def __init__(self, data, regions=None, title='', init_row=0, is_symmetric=False):
        BasePlotter1D.__init__(self, title=title)
        if regions is None:
            regions = []
            for i in range(len(data)):
                regions.append(GenomicRegion(chromosome='', start=i, end=i))

        self.init_row = init_row
        self.data = data
        self.sr = None
        self.da_sub = None
        self.regions = regions
        self.current_region = None
        self.line = None
        self.current_ix = init_row
        self.current_cutoff = None
        self.cutoff_line = None
        self.cutoff_line_mirror = None
        self.is_symmetric = is_symmetric

    def _new_region(self, region):
        self.current_region = region
        self.da_sub, self.sr = sub_data_regions(self.data, self.regions, region)

    def _plot(self, region=None, cax=None):
        self._new_region(region)
        bin_coords = [(x.start - 1) for x in self.sr]
        ds = self.da_sub[self.init_row]
        self.line, = self.ax.plot(bin_coords, ds)
        if not self.is_symmetric:
            self.current_cutoff = (self.ax.get_ylim()[1] - self.ax.get_ylim()[0]) / 2 + self.ax.get_ylim()[0]
        else:
            self.current_cutoff = self.ax.get_ylim()[1]/ 2
        self.ax.axhline(0.0, linestyle='dashed', color='grey')
        self.cutoff_line = self.ax.axhline(self.current_cutoff, color='r')
        if self.is_symmetric:
            self.cutoff_line_mirror = self.ax.axhline(-1*self.current_cutoff, color='r')
        self.ax.set_ylim((np.nanmin(ds), np.nanmax(ds)))

    def update(self, ix=None, cutoff=None, region=None, update_canvas=True):
        if region is not None:
            self._new_region(region)

        if ix is not None and ix != self.current_ix:
            ds = self.da_sub[ix]
            self.current_ix = ix
            self.line.set_ydata(ds)
            self.ax.set_ylim((np.nanmin(ds), np.nanmax(ds)))

            if cutoff is None:
                if not self.is_symmetric:
                    self.update(cutoff=(self.ax.get_ylim()[1]-self.ax.get_ylim()[0])/2 + self.ax.get_ylim()[0],
                                update_canvas=False)
                else:
                    self.update(cutoff=self.ax.get_ylim()[1] / 2, update_canvas=False)

            if update_canvas:
                self.fig.canvas.draw()

        if cutoff is not None and cutoff != self.current_cutoff:
            if self.is_symmetric:
                self.current_cutoff = abs(cutoff)
            else:
                self.current_cutoff = cutoff
            self.cutoff_line.set_ydata(self.current_cutoff)
            if self.is_symmetric:
                self.cutoff_line_mirror.set_ydata(-1*self.current_cutoff)

            if update_canvas:
                self.fig.canvas.draw()


class TADtoolPlot(object):
    def __init__(self, hic_matrix, regions=None, data=None, window_sizes=None, norm='lin', max_dist=3000000,
                 max_percentile=99.99, algorithm='insulation', matrix_colormap=None,
                 data_colormap=None, log_data=True):
        self.hic_matrix = hic_matrix
        if regions is None:
            regions = []
            for i in range(hic_matrix.shape[0]):
                regions.append(GenomicRegion(chromosome='', start=i, end=i))
        self.regions = regions
        self.norm = norm
        self.fig = None
        self.max_dist = max_dist
        self.algorithm = algorithm
        self.svmax = None
        self.min_value = np.nanmin(self.hic_matrix[np.nonzero(self.hic_matrix)])
        self.min_value_data = None
        self.hic_plot = None
        self.tad_plot = None
        self.data_plot = None
        self.line_plot = None
        self.sdata = None
        self.data_ax = None
        self.line_ax = None
        self.da = None
        self.ws = None
        self.current_window_size = None
        self.window_size_text = None
        self.tad_cutoff_text = None
        self.max_percentile = max_percentile
        self.tad_regions = None
        self.current_da_ix = None
        self.button_save_tads = None
        self.button_save_vector = None
        self.button_save_matrix = None

        self.log_data = log_data

        if algorithm == 'insulation':
            self.tad_algorithm = insulation_index
            self.tad_calling_algorithm = call_tads_insulation_index
            self.is_symmetric = False
            if matrix_colormap is None:
                self.matrix_colormap = LinearSegmentedColormap.from_list('myreds', ['white', 'red'])
            if data_colormap is None:
                self.data_plot_color = 'plasma'
        elif algorithm == 'ninsulation':
            self.tad_algorithm = normalised_insulation_index
            self.tad_calling_algorithm = call_tads_insulation_index
            self.is_symmetric = True
            if matrix_colormap is None:
                self.matrix_colormap = LinearSegmentedColormap.from_list('myreds', ['white', 'red'])
            if data_colormap is None:
                self.data_plot_color = LinearSegmentedColormap.from_list('myreds', ['blue', 'white', 'red'])
        elif algorithm == 'directionality':
            self.tad_algorithm = directionality_index
            self.tad_calling_algorithm = call_tads_directionality_index
            self.is_symmetric = True
            if matrix_colormap is None:
                self.matrix_colormap = LinearSegmentedColormap.from_list('myreds', ['white', 'red'])
            if data_colormap is None:
                self.data_plot_color = LinearSegmentedColormap.from_list('myreds', ['blue', 'white', 'red'])

        if data is None:
            self.da, self.ws = data_array(hic_matrix=self.hic_matrix, regions=self.regions,
                                          tad_method=self.tad_algorithm, window_sizes=window_sizes)
        else:
            self.da = data
            if window_sizes is None:
                raise ValueError("window_sizes parameter cannot be None when providing data!")
            self.ws = window_sizes

    def vmax_slider_update(self, val):
        self.hic_plot.set_clim(self.min_value, val)

    def data_slider_update(self, val):
        if self.is_symmetric:
            self.data_plot.set_clim(-1*val, val)
        else:
            self.data_plot.set_clim(self.min_value_data, val)

    def on_click_save_tads(self, event):
        tk.Tk().withdraw()  # Close the root window
        save_path = filedialog.asksaveasfilename()
        if save_path is not None:
            with open(save_path, 'w') as o:
                for region in self.tad_regions:
                    o.write("%s\t%d\t%d\n" % (region.chromosome, region.start-1, region.end))

    def on_click_save_vector(self, event):
        tk.Tk().withdraw()  # Close the root window
        save_path = filedialog.asksaveasfilename()
        if save_path is not None:
            da_sub = self.da[self.current_da_ix]
            with open(save_path, 'w') as o:
                for i, region in enumerate(self.regions):
                    o.write("%s\t%d\t%d\t.\t%e\n" % (region.chromosome, region.start-1, region.end, da_sub[i]))

    def on_click_save_matrix(self, event):
        tk.Tk().withdraw()  # Close the root window
        save_path = filedialog.asksaveasfilename()
        if save_path is not None:
            with open(save_path, 'w') as o:
                # write regions
                for i, region in enumerate(self.regions):
                    o.write("%s:%d-%d" % (region.chromosome, region.start-1, region.end))
                    if i < len(self.regions)-1:
                        o.write("\t")
                    else:
                        o.write("\n")

                # write matrix
                n_rows = self.da.shape[0]
                n_cols = self.da.shape[1]
                for i in range(n_rows):
                    window_size = self.ws[i]
                    o.write("%d\t" % window_size)

                    for j in range(n_cols):
                        o.write("%e" % self.da[i, j])
                        if j < n_cols-1:
                            o.write("\t")
                        else:
                            o.write("\n")

    def plot(self, region=None):
        # set up plotting grid
        self.fig = plt.figure(figsize=(10, 10))

        # main plots
        grid_size = (32, 15)
        hic_vmax_slider_ax = plt.subplot2grid(grid_size, (0, 0), colspan=13)
        hic_ax = plt.subplot2grid(grid_size, (1, 0), rowspan=9, colspan=13)
        hp_cax = plt.subplot2grid(grid_size, (1, 14), rowspan=9, colspan=1)
        tad_ax = plt.subplot2grid(grid_size, (10, 0), rowspan=1, colspan=13, sharex=hic_ax)
        line_ax = plt.subplot2grid(grid_size, (12, 0), rowspan=6, colspan=13, sharex=hic_ax)
        line_cax = plt.subplot2grid(grid_size, (12, 13), rowspan=6, colspan=2)
        data_vmax_slider_ax = plt.subplot2grid(grid_size, (19, 0), colspan=13)
        data_ax = plt.subplot2grid(grid_size, (20, 0), rowspan=9, colspan=13, sharex=hic_ax)
        da_cax = plt.subplot2grid(grid_size, (20, 14), rowspan=9, colspan=1)

        # buttons
        save_tads_ax = plt.subplot2grid(grid_size, (31, 0), rowspan=1, colspan=4)
        self.button_save_tads = Button(save_tads_ax, 'Save TADs')
        self.button_save_tads.on_clicked(self.on_click_save_tads)
        save_vector_ax = plt.subplot2grid(grid_size, (31, 5), rowspan=1, colspan=4)
        self.button_save_vector = Button(save_vector_ax, 'Save current values')
        self.button_save_vector.on_clicked(self.on_click_save_vector)
        save_matrix_ax = plt.subplot2grid(grid_size, (31, 10), rowspan=1, colspan=4)
        self.button_save_matrix = Button(save_matrix_ax, 'Save matrix')
        self.button_save_matrix.on_clicked(self.on_click_save_matrix)

        # add subplot content
        max_value = np.nanpercentile(self.hic_matrix, self.max_percentile)
        init_value = .2*max_value

        # HI-C VMAX SLIDER
        self.svmax = Slider(hic_vmax_slider_ax, 'vmax', self.min_value, max_value, valinit=init_value, color='grey')
        self.svmax.on_changed(self.vmax_slider_update)

        # HI-C
        self.hic_plot = HicPlot(self.hic_matrix, self.regions, max_dist=self.max_dist, norm=self.norm,
                                vmax=init_value, vmin=self.min_value, colormap=self.matrix_colormap)
        self.hic_plot.plot(region, ax=hic_ax, cax=hp_cax)

        # generate data array
        self.min_value_data = np.nanmin(self.da[np.nonzero(self.da)])
        max_value_data = np.nanpercentile(self.da, self.max_percentile)
        init_value_data = .5*max_value_data

        # LINE PLOT
        da_ix = int(self.da.shape[0]/2)
        self.current_da_ix = da_ix
        self.line_plot = DataLinePlot(self.da, regions=self.regions, init_row=da_ix, is_symmetric=self.is_symmetric)
        self.line_plot.plot(region, ax=line_ax)
        self.line_ax = line_ax

        # line info
        self.current_window_size = self.ws[da_ix]
        line_cax.text(.1, .8, 'Window size', fontweight='bold')
        self.window_size_text = line_cax.text(.3, .6, str(self.current_window_size))
        line_cax.text(.1, .4, 'TAD cutoff', fontweight='bold')
        self.tad_cutoff_text = line_cax.text(.3, .2, "%.5f" % self.line_plot.current_cutoff)
        line_cax.axis('off')

        # TAD PLOT
        self.tad_regions = self.tad_calling_algorithm(self.da[da_ix], self.line_plot.current_cutoff, self.regions)

        self.tad_plot = TADPlot(self.tad_regions)
        self.tad_plot.plot(region=region, ax=tad_ax)

        # DATA ARRAY
        self.data_plot = DataArrayPlot(self.da, self.ws, self.regions, vmax=init_value_data,
                                       colormap=self.data_plot_color, current_window_size=self.ws[da_ix],
                                       log_y=self.log_data)
        self.data_plot.plot(region, ax=data_ax, cax=da_cax)

        # DATA ARRAY SLIDER
        if self.is_symmetric:
            self.sdata = Slider(data_vmax_slider_ax, 'vmax', 0.0, max_value_data,
                                valinit=init_value_data, color='grey')
        else:
            self.sdata = Slider(data_vmax_slider_ax, 'vmax', self.min_value_data, max_value_data,
                                valinit=init_value_data, color='grey')

        self.sdata.on_changed(self.data_slider_update)
        self.data_slider_update(init_value_data)

        # clean up
        hic_ax.xaxis.set_visible(False)
        line_ax.xaxis.set_visible(False)

        # enable hover
        self.data_ax = data_ax
        cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

        return self.fig, (hic_vmax_slider_ax, hic_ax, line_ax, data_ax, hp_cax, da_cax)

    def on_click(self, event):
        if event.inaxes == self.data_ax or event.inaxes == self.line_ax:
            if event.inaxes == self.data_ax:
                ws_ix = bisect_left(self.ws, event.ydata)
                self.current_window_size = self.ws[ws_ix]
                self.current_da_ix = ws_ix

                self.data_plot.update(window_size=self.ws[ws_ix])
                self.line_plot.update(ix=ws_ix, update_canvas=False)
                self.tad_cutoff_text.set_text("%.5f" % self.line_plot.current_cutoff)
                self.window_size_text.set_text(str(self.current_window_size))
            elif event.inaxes == self.line_ax:
                if self.is_symmetric:
                    self.line_plot.update(cutoff=abs(event.ydata), update_canvas=False)
                else:
                    self.line_plot.update(cutoff=abs(event.ydata), update_canvas=False)
                self.tad_cutoff_text.set_text("%.5f" % self.line_plot.current_cutoff)

            # update TADs
            self.tad_regions = self.tad_calling_algorithm(self.da[self.current_da_ix], self.line_plot.current_cutoff,
                                                          self.regions)
            self.tad_plot.update(self.tad_regions)
            self.fig.canvas.draw()
