import tadtool.tad as tad
import tadtool.plot as tp

# regions = tad.HicRegionFileReader().regions("/Users/kkruse1/Desktop/kaic-demo/tadtool/regions_chr12_20-35Mb.txt")
# matrix = tad.HicMatrixFileReader().matrix("/Users/kkruse1/Desktop/kaic-demo/tadtool/chr12_20-35Mb.txt")


# da, ws = tad.data_array(matrix, regions)
# dap = tp.DataArrayPlot(da, ws, regions)
# fig, ax = dap.plot('chr12:25000000-30000000')
# fig.show()
regions = tad.HicRegionFileReader().regions("/Users/kkruse1/Desktop/kaic-demo/tadtool/chr12.10kb.regions.txt")
matrix = tad.HicMatrixFileReader().matrix("/Users/kkruse1/Desktop/kaic-demo/tadtool/chr12.10kb.txt")

# regions = tad.HicRegionFileReader().regions("/Users/kkruse1/Desktop/kaic-demo/tadtool/regions.txt")
# matrix = tad.HicMatrixFileReader().matrix("/Users/kkruse1/Desktop/kaic-demo/tadtool/chr12.mat")

ii = tp.insulation_index(matrix, regions, window_size=500000)
r = tad.call_tads_insulation_index(ii, 0.02, regions)

tad_plot = tp.TADtoolPlot(matrix, regions, norm='lin', max_dist=3000000, algorithm='insulation')
fig, axes = tad_plot.plot('chr12:31000000-34000000')
fig.show()

