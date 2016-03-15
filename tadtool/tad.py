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
