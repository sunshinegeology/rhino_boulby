import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pltpatches
import scipy.stats as scst


def goodness(sample, fit):
    return np.sqrt(np.mean([(sample[i]-fit[i])**2 for i in range(len(sample))]))
    

def histograms(fnames, prop_name):
    # csmp_fdata.txt: 'fracture' --> idx --> status --> ah --> rc --> nominal pressure --> normal gap --> shear gap --> rigid aperture --> mechanical aperture
    for fname in fnames:
        data = np.loadtxt(fname)
        
        # hist
        no_bins = 20
        hist, bin_edges = np.histogram(data, no_bins, normed=True)
        bin_centers = (bin_edges[:-1]+bin_edges[1:])/2
        
        # fit
        gfits = []
        dists = [scst.lognorm, scst.norm]
        for dist in dists:
            param = dist.fit(data)
            x_fitted = np.linspace(0., max(data), 100)
            gfits += [goodness(hist, dist.pdf(bin_centers, *param[:-2], loc=param[-2], scale=param[-1]))]
        dist = dists[np.argmin(gfits)]
        param = dist.fit(data)
        pdf_fitted = dist.pdf(x_fitted, *param[:-2], loc=param[-2], scale=param[-1])
        
        fig, ax = plt.subplots()
        color = 'blue'
        ax.plot(x_fitted, pdf_fitted, color=color, lw=1) 
        ax.plot(bin_centers, hist, marker='x', ls='None', color=color)
        #ax.hist(skewness, no_bins, normed=1)
        ax.set_xlabel(prop_name)
        ax.set_ylabel('Probability Density')
        plt.show()
        
if __name__ == "__main__":
    histograms(['skewness_east.txt', 'skewness_north.txt'], 'Skewness')