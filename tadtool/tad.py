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
