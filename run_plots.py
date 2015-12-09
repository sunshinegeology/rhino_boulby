import histograms

histograms.histograms(['skewness_east.txt', 'skewness_north.txt'], 'Skewness')
histograms.histograms(['aratio_east.txt', 'aratio_north.txt'], 'Aspect Ratio')

histograms.histograms(['skewness_east_observation.txt', 'skewness_north_observation.txt'], 'Skewness observation')