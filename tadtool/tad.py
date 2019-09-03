from __future__ import division
import numpy as np
import progressbar
import itertools
from future.utils import string_types
import warnings

import logging
logger = logging.getLogger(__name__)


class GenomicRegion(object):
    """
    Class representing a genomic region.

    .. attribute:: chromosome

        Name of the chromosome this region is located on

    .. attribute:: start

        Start position of the region in base pairs

    .. attribute:: end

        End position of the region in base pairs

    .. attribute:: strand

        Strand this region is on (+1, -1)

    .. attribute:: ix

        Index of the region in the context of all genomic
        regions.

    """

    def __init__(self, start, end, chromosome=None, ix=None):
        """
        Initialize this object.

        :param start: Start position of the region in base pairs
        :param end: End position of the region in base pairs
        :param chromosome: Name of the chromosome this region is located on
        :param ix: Index of the region in the context of all genomic
                   regions.
        """
        self.start = start
        self.end = end
        self.chromosome = chromosome
        self.ix = ix

    @classmethod
    def from_string(cls, region_string):
        """
        Convert a string into a :class:`~GenomicRegion`.

        This is a very useful convenience function to quickly
        define a :class:`~GenomicRegion` object from a descriptor
        string.

        :param region_string: A string of the form
                              <chromosome>[:<start>-<end>[:<strand>]]
                              (with square brackets indicating optional
                              parts of the string). If any optional
                              part of the string is omitted, intuitive
                              defaults will be chosen.
        :return: :class:`~GenomicRegion`
        """
        chromosome = None
        start = None
        end = None

        # strip whitespace
        no_space_region_string = "".join(region_string.split())
        fields = no_space_region_string.split(':')

        if len(fields) > 3:
            raise ValueError("Genomic range string must be of the form <chromosome>[:<start>-<end>:[<strand>]]")

        # there is chromosome information
        if len(fields) > 0:
            chromosome = fields[0]

        # there is range information
        if len(fields) > 1 and fields[1] != '':
            start_end_bp = fields[1].split('-')
            if len(start_end_bp) > 0:
                try:
                    start = int(start_end_bp[0])
                except ValueError:
                    raise ValueError("Start of genomic range must be integer")

            if len(start_end_bp) > 1:
                try:
                    end = int(start_end_bp[1])
                except ValueError:
                    raise ValueError("End of genomic range must be integer")

                if not end > start:
                    raise ValueError("The end coordinate must be bigger than the start.")

        return cls(start=start, end=end, chromosome=chromosome)

    def overlaps(self, region):
        """
        Check if this region overlaps with the specified region.

        :param region: :class:`~GenomicRegion` object or string
        """
        if isinstance(region, str):
            region = GenomicRegion.from_string(region)

        if region.chromosome != self.chromosome:
            return False

        if region.start <= self.end and region.end >= self.start:
            return True
        return False

    def contains(self, region):
        """
        Check if the specified region is completely contained in this region.

        :param region: :class:`~GenomicRegion` object or string
        """
        if isinstance(region, str):
            region = GenomicRegion.from_string(region)

        if region.chromosome != self.chromosome:
            return False

        if region.start >= self.start and region.end <= self.end:
            return True
        return False

    def _equals(self, region):
        if region.chromosome != self.chromosome:
            return False
        if region.start != self.start:
            return False
        if region.end != self.end:
            return False
        return True

    def __eq__(self, other):
        return self._equals(other)

    def __ne__(self, other):
        return not self._equals(other)

    def copy(self):
        """
        Make a copy of this GenomicRegion object.
        """
        return GenomicRegion(chromosome=self.chromosome, start=self.start, end=self.end)


def sub_regions(regions, region):
    """
    Get regions from a list the overlap with another region.

    :param regions: List of :class:`~GenomicRegion`
    :param region: :class:`~GenomicRegion` used for overlap calculation
    """
    if isinstance(region, string_types):
        region = GenomicRegion.from_string(region)

    sr = []
    start_ix = None
    end_ix = None
    for i, r in enumerate(regions):
        if r.chromosome == region.chromosome and r.start <= region.end and r.end >= region.start:
            if start_ix is None:
                start_ix = i
            end_ix = i
            sr.append(r.copy())
        else:
            if end_ix is not None:
                break

    if start_ix is None or end_ix is None:
        raise ValueError("Region not found in dataset! {}:{}-{}".format(region.chromosome, region.start, region.end))

    return sr, start_ix, end_ix


def sub_matrix_regions(hic_matrix, regions, region):
    """
    Get a square sub Hi-C matrix that overlaps a given region.

    :param hic_matrix: A square numpy array
    :param regions: List of :class:`~GenomicRegion`
    :param region: :class:`~GenomicRegion` used for overlap calculation
    """
    sr, start_ix, end_ix = sub_regions(regions, region)

    if start_ix is None:
        return np.empty((0, 0)), sr

    return np.copy(hic_matrix[start_ix:end_ix+1, start_ix:end_ix+1]), sr


def sub_data_regions(data, regions, region):
    """
    Get a sub data matrix that overlaps a given region.

    :param data: a numpy array
    :param regions: List of :class:`~GenomicRegion`
    :param region: :class:`~GenomicRegion` used for overlap calculation
    """
    sr, start_ix, end_ix = sub_regions(regions, region)

    if start_ix is None:
        return np.empty((0, 0)), sr

    if not isinstance(data, np.ndarray):
        data = np.array(data)

    return np.copy(data[:, start_ix:end_ix+1]), sr


def sub_vector_regions(data, regions, region):
    """
    Get a sub-vector that overlaps a given region

    :param data: a numpy vector
    :param regions: List of :class:`~GenomicRegion`
    :param region: :class:`~GenomicRegion` used for overlap calculation
    """
    sr, start_ix, end_ix = sub_regions(regions, region)

    if start_ix is None:
        return np.empty((0, 0)), sr

    if not isinstance(data, np.ndarray):
        data = np.array(data)

    return np.copy(data[start_ix:end_ix+1]), sr


def load_matrix(file_name, size=None, sep=None, square=True, ix_converter=None, region_range=None):
    try:  # numpy binary format
        m = np.load(file_name)
        if region_range is not None:
            s, e = region_range
            m = m[s:e+1, s:e+1].copy()
    except (IOError, ValueError):  # not an .npy file

        # check number of fields in file
        with open(file_name, 'r') as f:
            line = f.readline()
            while line.startswith('#'):
                line = f.readline()
            line = line.rstrip()
            n_fields = len(line.split(sep))

        if n_fields > 3 or not square:  # square matrix format
            m = np.loadtxt(file_name)
            if region_range is not None:
                s, e = region_range
                m = m[s:e + 1, s:e + 1].copy()
        else:
            if size is None:
                if region_range is not None:
                    size = region_range[1] - region_range[0] + 1
                else:
                    raise ValueError("Must provide matrix size when importing from sparse matrix notation "
                                     "(HicPro, etc)")

            m = np.zeros((size, size))
            with open(file_name, 'r') as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    line = line.rstrip()
                    fields = line.split(sep)
                    if ix_converter is None:
                        source, sink, weight = int(fields[0]), int(fields[1]), float(fields[2])
                    else:
                        source = ix_converter[fields[0]]
                        sink = ix_converter[fields[1]]
                        weight = float(fields[2])

                    if region_range is not None:
                        if source < region_range[0] or source > region_range[1]:
                            continue
                        if sink < region_range[0] or sink > region_range[1]:
                            continue
                        source = source - region_range[0]
                        sink = sink - region_range[0]

                    m[source, sink] = weight
                    m[sink, source] = weight

    if square and m.shape[0] != m.shape[1]:
        raise ValueError("Matrix dimensions do not match! ({})".format(m.shape))

    return m


def load_chromosome_matrix(file_name, regions, chromosome, **kwargs):
    region_range = [None, 0]
    for i, region in enumerate(regions):
        if region.chromosome == chromosome:
            if region_range[0] is None:
                region_range[0] = i
            region_range[0] = min(region_range[0], i)
            region_range[1] = max(region_range[1], i)
    logger.debug('Region range for {}: {} - {}'.format(chromosome, region_range[0], region_range[1]))
    return load_matrix(file_name, region_range=region_range, **kwargs)


def load_regions(file_name, sep=None):
    regions = []
    ix_converter = None
    with open(file_name, 'r') as f:
        for i, line in enumerate(f):
            line = line.rstrip()
            fields = line.split(sep)
            if len(fields) > 2:
                chromosome = fields[0]
                start = int(fields[1]) + 1
                end = int(fields[2])
                ix = i
                if len(fields) > 3 and fields[3] != '.':  # HicPro
                    if ix_converter is None:
                        ix_converter = dict()
                    if fields[3] in ix_converter:
                        raise ValueError("name column in region BED must "
                                         "only contain unique values! ({})".format(fields[3]))
                    ix_converter[fields[3]] = ix
                regions.append(GenomicRegion(chromosome=chromosome, start=start, end=end, ix=ix))
    return regions, ix_converter


def _get_boundary_distances(regions):
    """
    Return the distances of each region to the boundaries of its chromosome.
    """
    n_bins = len(regions)
    # find distances to chromosome boundaries in bins
    boundary_dist = np.zeros(n_bins, dtype=int)
    last_chromosome = None
    last_chromosome_index = 0
    for i, region in enumerate(regions):
        chromosome = region.chromosome
        if last_chromosome is not None and chromosome != last_chromosome:
            chromosome_length = i-last_chromosome_index
            for j in range(chromosome_length):
                boundary_dist[last_chromosome_index+j] = min(j, i-last_chromosome_index-1-j)
            last_chromosome_index = i
        last_chromosome = chromosome
    chromosome_length = n_bins-last_chromosome_index
    for j in range(chromosome_length):
        boundary_dist[last_chromosome_index+j] = min(j, n_bins-last_chromosome_index-1-j)

    return boundary_dist


def directionality_index(hic, regions=None, window_size=2000000):
    """
    Calculate the directionality index as in Dixon et al. (2012).

    :param hic: A Hi-C matrix in the form of a square numpy array
    :param regions: A list of :class:`~GenomicRegion`s - if omitted, will create a dummy list
    :param window_size: A window size in base pairs
    """
    if regions is None:
        for i in range(hic.shape[0]):
            regions.append(GenomicRegion(chromosome='', start=i, end=i))

    bin_size = regions[0].end-regions[0].start+1
    bin_distance = int(window_size/bin_size)
    if window_size % bin_size > 0:
        bin_distance += 1

    n_bins = len(regions)
    boundary_dist = _get_boundary_distances(regions)

    left_sums = np.zeros(n_bins)
    right_sums = np.zeros(n_bins)
    for source in range(hic.shape[0]):
        for sink in range(source, hic.shape[1]):
            try:
                if hic.mask[source, sink]:
                    continue
            except AttributeError:
                pass

            weight = hic[source, sink]

            if source == sink:
                continue
            if sink - source <= bin_distance:
                if boundary_dist[sink] >= sink-source:
                    left_sums[sink] += weight
                if boundary_dist[source] >= sink-source:
                    right_sums[source] += weight

    di = np.zeros(n_bins)
    for i in range(n_bins):
        A = left_sums[i]
        B = right_sums[i]
        E = (A+B)/2
        if E != 0 and B-A != 0:
            di[i] = ((B-A)/abs(B-A)) * ((((A-E)**2)/E) + (((B-E)**2)/E))

    return di


def masked_matrix(matrix, all_zero=False):
    """
    Returns masked version of HicMatrix. By default, all entries in zero-count
    rows and columns are masked.

    :param matrix: A numpy 2D matrix
    :param all_zero: Mask ALL zero-count entries
    :returns: MaskedArray with zero entries masked
    """
    if all_zero:
        return np.ma.MaskedArray(matrix, mask=np.isclose(matrix, 0.))
    col_zero = np.isclose(np.sum(matrix, axis=0), 0.)
    row_zero = np.isclose(np.sum(matrix, axis=1), 0.)
    mask = np.zeros(matrix.shape, dtype=np.bool_)
    mask[:, col_zero] = np.True_
    mask[row_zero, :] = np.True_
    return np.ma.MaskedArray(matrix, mask=mask)


def kth_diag_indices(n, k):
    """
    Return indices of bins k steps away from the diagonal.
    (from http://stackoverflow.com/questions/10925671/numpy-k-th-diagonal-indices)
    """

    rows, cols = np.diag_indices(n)
    if k < 0:
        return rows[:k], cols[-k:]
    elif k > 0:
        return rows[k:], cols[:-k]
    else:
        return rows, cols


def impute_missing_bins(hic_matrix, regions=None, per_chromosome=True, stat=np.ma.mean):
    """
    Impute missing contacts in a Hi-C matrix.

    For inter-chromosomal data uses the mean of all inter-chromosomal contacts,
    for intra-chromosomal data uses the mean of intra-chromosomal counts at the corresponding diagonal.

    :param hic_matrix: A square numpy array
    :param regions: A list of :class:`~GenomicRegion`s - if omitted, will create a dummy list
    :param per_chromosome: Do imputation on a per-chromosome basis (recommended)
    :param stat: The aggregation statistic to be used for imputation, defaults to the mean.
    """
    if regions is None:
        for i in range(hic_matrix.shape[0]):
            regions.append(GenomicRegion(chromosome='', start=i, end=i))

    chr_bins = dict()
    for i, region in enumerate(regions):
        if region.chromosome not in chr_bins:
            chr_bins[region.chromosome] = [i, i]
        else:
            chr_bins[region.chromosome][1] = i

    n = len(regions)
    if not hasattr(hic_matrix, "mask"):
        hic_matrix = masked_matrix(hic_matrix)

    imputed = hic_matrix.copy()
    if per_chromosome:
        for c_start, c_end in chr_bins.itervalues():
            # Correcting intrachromoc_startmal contacts by mean contact count at each diagonal
            for i in range(c_end - c_start):
                ind = kth_diag_indices(c_end - c_start, -i)
                diag = imputed[c_start:c_end, c_start:c_end][ind]
                diag[diag.mask] = stat(diag)
                imputed[c_start:c_end, c_start:c_end][ind] = diag
            # Correcting interchromoc_startmal contacts by mean of all contact counts between
            # each set of chromoc_startmes
            for other_start, other_end in chr_bins.itervalues():
                # Only correct upper triangle
                if other_start <= c_start:
                    continue
                inter = imputed[c_start:c_end, other_start:other_end]
                inter[inter.mask] = stat(inter)
                imputed[c_start:c_end, other_start:other_end] = inter
    else:
        for i in range(n):
            diag = imputed[kth_diag_indices(n, -i)]
            diag[diag.mask] = stat(diag)
            imputed[kth_diag_indices(n, -i)] = diag
    # Copying upper triangle to lower triangle
    imputed[np.tril_indices(n)] = imputed.T[np.tril_indices(n)]
    return imputed


def _apply_sliding_func(a, window, func=np.ma.mean):
    """
    Apply function on a sliding window over an array, ignoring Numpy NaN values.

    :param a: Numpy array on which function is applied
    :param window: The sliding window is i - window:i + window + 1
                   so total window is twice this parameter.
    :param func: Aggregation function to apply
    """
    out = np.empty(len(a))
    for i in range(len(a)):
        window_start = max(0, i - window)
        window_end = min(len(a), i + window + 1)
        cur_window = np.array(a[window_start:window_end])
        out[i] = func(cur_window[~np.isnan(cur_window)])
    return out


def insulation_index(hic_matrix, regions=None, window_size=500000, relative=False, mask_thresh=.5,
                     aggr_func=np.nanmean, impute_missing=False, normalize=False,
                     normalization_window=None, gradient=False):
    """
    Calculate the insulation index as in Crane et al. (2015)

    :param hic_matrix: A square numpy array
    :param regions: A list of :class:`~GenomicRegion`s - if omitted, will create a dummy list
    :param window_size: A window size in base pairs
    :param relative: Calculate insulation index relative to upstream and downstream contact strength
    :param mask_thresh: Set the threshold for masked value threshold
    :param aggr_func: Aggregation function used for index calculation, defaults to the mean
    :param impute_missing: Impute missing values, defaults to False
    :param normalize: Normalize index at each region to the mean number of contacts across a larger region
                      and log2-transform data
    :param normalization_window: Window (in number of regions) to normalize contacts to (see normalize).
                                 If None, the whole chromosome will be used for normalization.
    :param gradient: Return the first derivative of the insulation index
    """
    if regions is None:
        for i in range(hic_matrix.shape[0]):
            regions.append(GenomicRegion(chromosome='', start=i, end=i))

    bin_size = regions[0].end-regions[0].start+1
    d = int(window_size/bin_size)
    if window_size % bin_size > 0:
        d += 1

    chr_bins = dict()
    for i, region in enumerate(regions):
        if region.chromosome not in chr_bins:
            chr_bins[region.chromosome] = [i, i]
        else:
            chr_bins[region.chromosome][1] = i

    if impute_missing:
        hic_matrix = impute_missing_bins(hic_matrix, regions=regions)

    is_masked = False
    if hasattr(hic_matrix, 'mask'):
        is_masked = True

    skipped = 0
    last_chromosome = None
    ins_by_chromosome = []
    for i, r in enumerate(regions):
        if r.chromosome != last_chromosome:
            last_chromosome = r.chromosome
            ins_by_chromosome.append(list())

        # too close to chromosome edge
        if (i - chr_bins[r.chromosome][0] < d or
                chr_bins[r.chromosome][1] - i <= d + 1):
            ins_by_chromosome[-1].append(np.nan)
            # ins_matrix[i] = np.nan
            continue

        # masked region
        if is_masked and hic_matrix.mask[i, i]:
            ins_by_chromosome[-1].append(np.nan)
            # ins_matrix[i] = np.nan
            continue

        up_rel_slice = (slice(i - d, i), slice(i - d, i))
        down_rel_slice = (slice(i + 1, i + d + 1), slice(i + 1, i + d + 1))
        ins_slice = (slice(i + 1, i + d + 1), slice(i - d, i))

        if (is_masked and ((relative and np.sum(hic_matrix.mask[up_rel_slice]) > d*d*mask_thresh) or
                (relative and np.sum(hic_matrix.mask[down_rel_slice]) > d*d*mask_thresh) or
                np.sum(hic_matrix.mask[ins_slice]) > d*d*mask_thresh)):
            # If too close to the edge of chromosome or
            # if more than half of the entries in this quadrant are masked (unmappable)
            # exclude it from the analysis
            skipped += 1
            ins_by_chromosome[-1].append(np.nan)
            # ins_matrix[i] = np.nan
            continue

        if not relative:
            if impute_missing and is_masked:
                s = aggr_func(hic_matrix[ins_slice].data)
                ins_by_chromosome[-1].append(s)
                # ins_matrix[i] = s
            else:
                s = aggr_func(hic_matrix[ins_slice])
                ins_by_chromosome[-1].append(s)
                # ins_matrix[i] = s
        else:
            if not impute_missing and not is_masked:
                s = (aggr_func(hic_matrix[ins_slice]) /
                     aggr_func(np.ma.dstack((hic_matrix[up_rel_slice],
                                             hic_matrix[down_rel_slice]))))
                # ins_matrix[i] = s
                ins_by_chromosome[-1].append(s)
            else:
                s = (aggr_func(hic_matrix[ins_slice].data) /
                     aggr_func(np.ma.dstack((hic_matrix[up_rel_slice].data,
                                             hic_matrix[down_rel_slice].data))))
                # ins_matrix[i] = s
                ins_by_chromosome[-1].append(s)

    if normalize:
        for i in range(len(ins_by_chromosome)):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if normalization_window is not None:
                    ins_by_chromosome[i] = np.ma.log2(ins_by_chromosome[i]/_apply_sliding_func(
                        ins_by_chromosome[i], normalization_window, func=np.nanmean))
                else:
                    ins_by_chromosome[i] = np.ma.log2(ins_by_chromosome[i]/np.nanmean(ins_by_chromosome[i]))
    if gradient:
        for i in range(len(ins_by_chromosome)):
            ins_by_chromosome[i] = np.gradient(ins_by_chromosome[i])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ins_matrix = np.array(list(itertools.chain.from_iterable(ins_by_chromosome)))

    return ins_matrix


def normalised_insulation_index(hic_matrix, regions=None, window_size=500000, normalisation_window=None):
    return insulation_index(hic_matrix, regions=regions, window_size=window_size, normalize=True,
                            normalization_window=normalisation_window)


def data_array(hic_matrix, regions, tad_method=insulation_index,
               window_sizes=None, **kwargs):
    """
    Calculate an index for a range of window sizes.

    :param hic_matrix: A square numpy array
    :param regions: A list of :class:`~GenomicRegion`s - if omitted, will create a dummy list
    :param tad_method: The method used for index calculation
    :param window_sizes: A list of window sizes in base pairs
    :param kwargs: Parameters passed to tad_method
    """
    if window_sizes is None:
        window_sizes = [int(i) for i in np.logspace(4, 6, 100)]

    if regions is None:
        for i in range(hic_matrix.shape[0]):
            regions.append(GenomicRegion(chromosome='', start=i, end=i))

    try:
        l = len(window_sizes)
    except TypeError:
        l = None

    pb = progressbar.ProgressBar(max_value=l)
    if l is not None:
        pb.start()

    ws = []
    da = []
    for i, window_size in enumerate(window_sizes):
        if l is not None:
            pb.update(i)
        ws.append(window_size)
        da.append(tad_method(hic_matrix, regions, window_size, **kwargs))

    if l is not None:
        pb.finish()

    return np.array(da), ws


def _border_type(i, values):
    """
    Returns border type:

    :param i: index of potential border
    :param values: insulation index values
    :return: if border: 1 or -1, else 0. 1 if derivative at border is >0, -1 otherwise
    """

    if i == 0:
        return 1
    if i == len(values)-1:
        return -1

    if values[i-1] <= 0 <= values[i+1]:
        return 1
    if values[i+1] <= 0 <= values[i-1]:
        return -1
    return 0


def call_tad_borders(ii_results, cutoff=0):
    """
    Calls TAD borders using the first derivative of the insulation index.

    :param ii_results: (raw) insulation index results, numpy vector
    :param cutoff: raw insulation index threshold for "true" TAD boundaries
    """
    tad_borders = []
    g = np.gradient(ii_results)
    for i in range(len(ii_results)):
        border_type = _border_type(i, g)
        if border_type == 1 and ii_results[i] <= cutoff:
            tad_borders.append(i)
    return tad_borders


def call_tads_insulation_index(ii_results, cutoff, regions=None):
    """
    Call TADs from insulation index vector.

    :param ii_results: (raw) insulation index results, numpy vector
    :param cutoff: raw insulation index threshold for "true" TADs
    :param regions: A list of :class:`~GenomicRegion`s - if omitted, will create a dummy list
    """
    if regions is None:
        regions = []
        for i in range(len(ii_results)):
            regions.append(GenomicRegion(chromosome='', start=i, end=i))

    if len(regions) != len(ii_results):
        raise ValueError("Insulation index results and regions must be the same size!")

    borders = call_tad_borders(ii_results, cutoff=cutoff)

    tad_regions = []
    previous = None
    for border in borders:
        if previous is not None:
            found_max = False
            for j in range(previous, border+1):
                if ii_results[j] >= cutoff:
                    found_max = True

            if found_max:
                first_region = regions[previous]
                second_region = regions[border]
                if first_region.chromosome != second_region.chromosome:
                    previous = border
                    continue
                tad_regions.append(GenomicRegion(chromosome=first_region.chromosome,
                                                 start=first_region.start, end=second_region.end))

        previous = border
    return tad_regions


def call_tads_directionality_index(di_results, cutoff, regions=None):
    """
    Call TADs from directionality index vector.

    :param di_results: (raw) directionality index results, numpy vector
    :param cutoff: raw directionality index threshold for "true" biases
    :param regions: A list of :class:`~GenomicRegion`s - if omitted, will create a dummy list
    """
    if regions is None:
        regions = []
        for i in range(len(di_results)):
            regions.append(GenomicRegion(chromosome='', start=i, end=i))

    tad_regions = []
    last_upstream_bias_region = None
    last_downstream_bias_region = None
    for i, value in enumerate(di_results):
        current_region = regions[i]
        if last_upstream_bias_region is not None and current_region.chromosome != last_upstream_bias_region.chromosome:
            tad_regions.append(GenomicRegion(chromosome=last_upstream_bias_region.chromosome,
                                             start=last_upstream_bias_region.start,
                                             end=regions[i-1].end))
            last_upstream_bias_region = None
            last_downstream_bias_region = None

        # look for upstream bias site
        if value >= cutoff:
            # side-by-side TAD
            if last_downstream_bias_region is not None and last_upstream_bias_region is not None:
                tad = GenomicRegion(chromosome=last_downstream_bias_region.chromosome,
                                    start=last_upstream_bias_region.start,
                                    end=last_downstream_bias_region.end)
                tad_regions.append(tad)
                last_upstream_bias_region = None
                last_downstream_bias_region = None

            # beginning of TAD
            if last_upstream_bias_region is None:
                last_upstream_bias_region = current_region

        if value <= -1*cutoff:
            # potential end of TAD
            if last_upstream_bias_region is not None:
                last_downstream_bias_region = current_region

        if abs(value) < cutoff:
            # back to zero
            if last_downstream_bias_region is not None and last_upstream_bias_region is not None:
                tad = GenomicRegion(chromosome=last_downstream_bias_region.chromosome,
                                    start=last_upstream_bias_region.start,
                                    end=last_downstream_bias_region.end)
                tad_regions.append(tad)
                last_upstream_bias_region = None
                last_downstream_bias_region = None
    return tad_regions
