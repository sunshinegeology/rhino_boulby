import numpy as np
import matplotlib.pyplot as plt


def two_histograms(fname_sim, fname_obs, ax):
    simulation_data = np.loadtxt(fname_sim)
    observation_data = np.loadtxt(fname_obs)
    bin_observation = 20
    bin_simulation = 20
    ax.hist(simulation_data, color = 'b', bins = bin_simulation, alpha = 1.0, normed = True, label = 'Simulation')
    ax.hist(observation_data, color = 'k', bins = bin_observation, alpha = 0.5, normed = True,  label = 'Observation')
    ax.legend(fontsize=12)
    


if __name__ == '__main__':
    # east
    # 2x2 plots
    f, axarr = plt.subplots(2, 2, figsize=(14,8))
    # width
    axarr[0, 0].set_title('Width')
    fname_sim = 'width_east.txt'
    fname_obs = 'width_east_observation.txt'
    two_histograms(fname_sim, fname_obs, axarr[0, 0])
    # height
    axarr[0, 1].set_title('Height')
    fname_sim = 'height_east.txt'
    fname_obs = 'height_east_observation.txt'
    two_histograms(fname_sim, fname_obs, axarr[0, 1])
    # aspect ratio
    axarr[1, 1].set_title('Aspect Ratio')
    fname_sim = 'aratio_east.txt'
    fname_obs = 'aratio_east_observation.txt'
    two_histograms(fname_sim, fname_obs, axarr[1, 1])
    # skewness ratio
    axarr[1, 0].set_title('Skewness')
    fname_sim = 'aratio_east.txt'
    fname_obs = 'aratio_east_observation.txt'
    two_histograms(fname_sim, fname_obs, axarr[1, 0])
    plt.savefig('histograms_comparison.png', dpi=300)
    plt.show()