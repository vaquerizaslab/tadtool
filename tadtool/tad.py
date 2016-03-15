from __future__ import division
import numpy as np


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
                    regions.append([chromosome, start, end])
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
            chromosome = region[0]
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
            regions.append(['', i, i])

    bin_size = regions[0][2]-regions[0][1]+1
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
            regions.append(['', i, i])

    chr_bins = dict()
    for i, region in enumerate(regions):
        if region[0] not in chr_bins:
            chr_bins[region[0]] = [i, i]
        else:
            chr_bins[region[0]][1] = i

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