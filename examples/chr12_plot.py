import tadtool.tad as tad
import tadtool.plot as tp

# load regions data set
regions = tad.HicRegionFileReader().regions("chr12_20-35Mb_regions.bed")

# load matrix
matrix = tad.HicMatrixFileReader().matrix("chr12_20-35Mb.matrix.txt")

# prepare plot
tad_plot = tp.TADtoolPlot(matrix, regions, norm='lin', max_dist=1000000, algorithm='insulation')
fig, axes = tad_plot.plot('chr12:31000000-34000000')

# show plot
import matplotlib.pyplot as plt
plt.show()

