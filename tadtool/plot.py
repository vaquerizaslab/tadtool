from __future__ import division, print_function
from abc import ABCMeta, abstractmethod
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator, Formatter, Locator
from matplotlib.widgets import Slider
import matplotlib.pyplot as plt
from tadtool.tad import GenomicRegion, sub_matrix_regions
import math
import copy
import numpy as np


class BasePlotter(object):

    __metaclass__ = ABCMeta

    def __init__(self, title):
        self._ax = None
        self.cax = None
        self.title = title

    @abstractmethod
    def _plot(self, region=None):
        raise NotImplementedError("Subclasses need to override _plot function")

    @abstractmethod
    def _refresh(self, region=None):
        raise NotImplementedError("Subclasses need to override _refresh function")

    @abstractmethod
    def plot(self, region=None):
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

    def plot(self, region=None, ax=None):
        if isinstance(region, basestring):
            region = GenomicRegion.from_string(region)
        if ax:
            self.ax = ax
        # set genome tick formatter
        self.ax.xaxis.set_major_formatter(GenomeCoordFormatter(region))
        self.ax.xaxis.set_major_locator(GenomeCoordLocator(nbins=5))
        self.ax.xaxis.set_minor_locator(MinorGenomeCoordLocator(n=5))
        self.ax.set_title(self.title)
        self._plot(region)
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
            for i in xrange(hic_matrix.shape[0]):
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

    def _plot(self, region=None):
        if region is None:
            raise ValueError("Cannot plot triangle plot for whole genome.")

        hm, sr = sub_matrix_regions(self.hic_matrix, self.regions, region)
        hm[np.tril_indices(hm.shape[0])] = np.nan
        # Remove part of matrix further away than max_dist
        if self.max_dist is not None:
            for i in xrange(hm.shape[0]):
                i_region = sr[i]
                for j in xrange(hm.shape[1]):
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
        self.ax.set_aspect(aspect="equal")
        self.ax.set_ylim(0, self.max_dist/2 if self.max_dist else 0.5*(region.end-region.start))
        # remove y ticks
        self.ax.set_yticks([])
        # hide background patch
        self.ax.patch.set_visible(False)
        if self.show_colorbar:
            self.add_colorbar(None)

    def set_clim(self, vmin, vmax):
        self.hicmesh.set_clim(vmin=vmin, vmax=vmax)
        if self.colorbar is not None:
            self.colorbar.set_clim(vmin=vmin, vmax=vmax)
            self.colorbar.draw_all()

    def _refresh(self, region=None):
        pass


class TADtoolPlot(object):
    def __init__(self, hic_matrix, regions=None, norm='lin', max_dist=3000000,
                 max_percentile=99.99):
        self.hic_matrix = hic_matrix
        if regions is None:
            regions = []
            for i in xrange(hic_matrix.shape[0]):
                regions.append(GenomicRegion(chromosome='', start=i, end=i))
        self.regions = regions
        self.norm = norm
        self.max_dist = max_dist
        self.svmax = None
        self.min_value = np.nanmin(self.hic_matrix[np.nonzero(self.hic_matrix)])
        self.hic_plot = None
        self.max_percentile = max_percentile

    def slider_update(self, val):
        self.hic_plot.set_clim(self.min_value, val)

    def plot(self, region=None):
        # set up plotting grid
        fig = plt.figure(figsize=(20, 4))
        ax1 = plt.subplot2grid((10, 1), (0, 0))
        ax2 = plt.subplot2grid((10, 1), (1, 0), rowspan=9)

        # add subplot content
        max_value = np.nanpercentile(self.hic_matrix, self.max_percentile)
        init_value = .7*max_value

        self.hic_plot = HicPlot(self.hic_matrix, self.regions, max_dist=self.max_dist, norm=self.norm,
                                vmax=init_value, vmin=self.min_value)
        self.svmax = Slider(ax1, 'vmax', self.min_value, max_value, valinit=init_value)
        self.svmax.on_changed(self.slider_update)

        self.hic_plot.plot(region, ax=ax2)

        return fig, (ax1, ax2)
