import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

def two_histograms(fname_sim, fname_obs, ax):
    simulation_data = np.loadtxt(fname_sim)
    observation_data = np.loadtxt(fname_obs)
    bin_observation = 20
    bin_simulation = 20
    ax.hist(simulation_data, color = 'b', bins = bin_simulation, alpha = 1.0, normed = True, label = 'Simulation')
    ax.hist(observation_data, color = 'k', bins = bin_observation, alpha = 0.5, normed = True,  label = 'Observation')
    ax.legend(fontsize=12)
    ax.text(0.5, 0.95, '{0} = {1:d}\n{2} = {3:d}'.format('$\mathit{N_S}$', len(simulation_data), '$\mathit{N_O}$', len(observation_data)), transform=ax.transAxes, va='top')


if __name__ == '__main__':
    # TODO north data too
    for ort in ['east', 'north']:
        # east
        # 2x2 plots
        f, axarr = plt.subplots(2, 2, figsize=(14,8))
        # width
        # axarr[0, 0].set_title('Width')
        axarr[0, 0].set_xlabel('Width (m)', fontsize = 14)
        axarr[0, 0].set_ylabel('Relative frequency', fontsize = 14)
        fname_sim = 'width_'+ort+'.txt'
        fname_obs = 'width_'+ort+'_observation.txt'
        two_histograms(fname_sim, fname_obs, axarr[0, 0])
        # height
        # axarr[0, 1].set_title('Height')
        axarr[0, 1].set_xlabel('Height (m)', fontsize = 14)
        fname_sim = 'height_'+ort+'.txt'
        fname_obs = 'height_'+ort+'_observation.txt'
        two_histograms(fname_sim, fname_obs, axarr[0, 1])
        # aspect ratio
        # axarr[1, 1].set_title('Aspect Ratio')
        axarr[1, 1].set_xlabel('Aspect Ratio', fontsize = 14)
        fname_sim = 'aratio_'+ort+'.txt'
        fname_obs = 'aratio_'+ort+'_observation.txt'
        two_histograms(fname_sim, fname_obs, axarr[1, 1])
        # skewness ratio
        # axarr[1, 0].set_title('Skewness')
        axarr[1, 0].set_xlabel('Skewness', fontsize = 14)
        axarr[1, 0].set_ylabel('Relative frequency', fontsize = 14)
        fname_sim = 'skewness_'+ort+'.txt'
        fname_obs = 'skewness_'+ort+'_observation.txt'
        two_histograms(fname_sim, fname_obs, axarr[1, 0])
        plt.savefig('histograms_comparison_'+ort+'.png', dpi=300)
        plt.show()