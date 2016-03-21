from __future__ import division
import numpy as np
import progressbar


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
        return GenomicRegion(chromosome=self.chromosome, start=self.start, end=self.end)


def sub_regions(regions, region):
    if isinstance(region, basestring):
        region = GenomicRegion.from_string(region)

    sr = []
    start_ix = None
    end_ix = None
    for i, r in enumerate(regions):
        if r.start <= region.end and r.end >= region.start:
            if start_ix is None:
                start_ix = i
            end_ix = i
            sr.append(r.copy())
        else:
            if end_ix is not None:
                break

    return sr, start_ix, end_ix


def sub_matrix_regions(hic_matrix, regions, region):
    sr, start_ix, end_ix = sub_regions(regions, region)

    if start_ix is None:
        return np.empty((0, 0)), sr

    return np.copy(hic_matrix[start_ix:end_ix+1, start_ix:end_ix+1]), sr


def sub_data_regions(data, regions, region):
    sr, start_ix, end_ix = sub_regions(regions, region)

    if start_ix is None:
        return np.empty((0, 0)), sr

    if not isinstance(data, np.ndarray):
        data = np.array(data)

    return np.copy(data[:, start_ix:end_ix+1]), sr


def sub_vector_regions(data, regions, region):
    sr, start_ix, end_ix = sub_regions(regions, region)

    if start_ix is None:
        return np.empty((0, 0)), sr

    if not isinstance(data, np.ndarray):
        data = np.array(data)

    return np.copy(data[start_ix:end_ix+1]), sr


class HicMatrixFileReader(object):
    def __init__(self, file_name=None):
        self.m = None

        if file_name is not None:
            self.load(file_name)

    def load(self, file_name):
        self.m = np.loadtxt(file_name)

    def matrix(self, file_name=None):
        if self.m is None and file_name is not None:
            self.load(file_name)
        else:
            raise ValueError("Either provide a file name or explicitly use 'load' before requesting matrix")
        return self.m


class HicRegionFileReader(object):
    def __init__(self, file_name=None, _separator="\t"):
        self.sep = _separator
        self.r = None
        if file_name is not None:
            self.load(file_name)

    def load(self, file_name):
        regions = []
        with open(file_name, 'r') as f:
            for line in f:
                line = line.rstrip()
                fields = line.split(self.sep)
                if len(fields) > 2:
                    chromosome = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    regions.append(GenomicRegion(chromosome=chromosome, start=start, end=end))
        self.r = regions

    def regions(self, file_name=None):
        if self.r is None and file_name is not None:
            self.load(file_name)
        else:
            raise ValueError("Either provide a file name or explicitly use 'load' before requesting regions")
        return self.r


def _get_boundary_distances(regions):
        n_bins = len(regions)
        # find distances to chromosome boundaries in bins
        boundary_dist = np.zeros(n_bins, dtype=int)
        last_chromosome = None
        last_chromosome_index = 0
        for i, region in enumerate(regions):
            chromosome = region.chromosome
            if last_chromosome is not None and chromosome != last_chromosome:
                chromosome_length = i-last_chromosome_index
                for j in xrange(chromosome_length):
                    boundary_dist[last_chromosome_index+j] = min(j, i-last_chromosome_index-1-j)
                last_chromosome_index = i
            last_chromosome = chromosome
        chromosome_length = n_bins-last_chromosome_index
        for j in xrange(chromosome_length):
            boundary_dist[last_chromosome_index+j] = min(j, n_bins-last_chromosome_index-1-j)

        return boundary_dist


def directionality_index(hic, regions=None, window_size=2000000):
    if regions is None:
        for i in xrange(hic.shape[0]):
            regions.append(GenomicRegion(chromosome='', start=i, end=i))

    bin_size = regions[0].end-regions[0].start+1
    bin_distance = int(window_size/bin_size)
    if window_size % bin_size > 0:
        bin_distance += 1

    n_bins = len(regions)
    boundary_dist = _get_boundary_distances(regions)

    left_sums = np.zeros(n_bins)
    right_sums = np.zeros(n_bins)
    for source in xrange(hic.shape[0]):
        for sink in xrange(source, hic.shape[1]):
            weight = hic[source, sink]

            if source == sink:
                continue
            if sink - source <= bin_distance:
                if boundary_dist[sink] >= sink-source:
                    left_sums[sink] += weight
                if boundary_dist[source] >= sink-source:
                    right_sums[source] += weight

    di = np.zeros(n_bins)
    for i in xrange(n_bins):
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
    # http://stackoverflow.com/questions/10925671/numpy-k-th-diagonal-indices
    rows, cols = np.diag_indices(n)
    if k < 0:
        return rows[:k], cols[-k:]
    elif k > 0:
        return rows[k:], cols[:-k]
    else:
        return rows, cols


def impute_missing_bins(hic_matrix, regions=None, per_chromosome=True, stat=np.ma.mean):
    if regions is None:
        for i in xrange(hic_matrix.shape[0]):
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


def insulation_index(hic_matrix, regions=None, window_size=2000000, relative=False, mask_thresh=.5,
                     aggr_func=np.nanmean, impute_missing=False):
    if regions is None:
        for i in xrange(hic_matrix.shape[0]):
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

    n = len(regions)
    ins_matrix = np.empty(n)
    skipped = 0
    for i, r in enumerate(regions):
        # too close to chromosome edge
        if (i - chr_bins[r.chromosome][0] < d or
                chr_bins[r.chromosome][1] - i <= d + 1):
            ins_matrix[i] = np.nan
            continue

        # masked region
        if is_masked and hic_matrix.mask[i, i]:
            ins_matrix[i] = np.nan
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
            ins_matrix[i] = np.nan
            continue

        if not relative:
            if impute_missing and is_masked:
                s = aggr_func(hic_matrix[ins_slice].data)
                ins_matrix[i] = s
            else:
                s = aggr_func(hic_matrix[ins_slice])
                ins_matrix[i] = s
        else:
            if not impute_missing and not is_masked:
                ins_matrix[i] = (aggr_func(hic_matrix[ins_slice]) /
                                 aggr_func(np.ma.dstack((hic_matrix[up_rel_slice],
                                                         hic_matrix[down_rel_slice]))))
            else:
                ins_matrix[i] = (aggr_func(hic_matrix[ins_slice].data) /
                                 aggr_func(np.ma.dstack((hic_matrix[up_rel_slice].data,
                                                         hic_matrix[down_rel_slice].data))))

    return ins_matrix


def data_array(hic_matrix, regions, tad_method=directionality_index,
               window_sizes=list(xrange(20000, 1000000, 20000)), **kwargs):
    if regions is None:
        for i in xrange(hic_matrix.shape[0]):
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


def call_tads_insulation_index(ii_results, cutoff, regions=None):
    if regions is None:
        for i in xrange(len(ii_results)):
            regions.append(GenomicRegion(chromosome='', start=i, end=i))

    if len(regions) != len(ii_results):
        raise ValueError("Insulation index results and regions must be the same size!")

    tad_regions = []
    current_tad_start = None
    current_tad_chromosome = None
    for i, value in enumerate(ii_results):
        current_region = regions[i]
        if current_tad_chromosome is not None and current_region.chromosome != current_tad_chromosome:
            tad_regions.append(GenomicRegion(chromosome=current_tad_chromosome, start=current_tad_start,
                                             end=regions[i-1].end))
            current_tad_chromosome = None
            current_tad_start = None

        if value >= cutoff:
            if current_tad_start is None:
                current_tad_start = current_region.start
                current_tad_chromosome = current_region.chromosome
        elif current_tad_start is not None:
                tad_regions.append(GenomicRegion(chromosome=current_tad_chromosome, start=current_tad_start,
                                                 end=regions[i-1].end))
                current_tad_chromosome = None
                current_tad_start = None

    return tad_regions
