import numpy as np
import matplotlib.pyplot as plt

simulation_data = np.loadtxt('width_east.txt')
observation_data = np.loadtxt('width_observation.txt')
bin_observation = 20
bin_simulation = 20
plt.hist(simulation_data, color = 'b', bins = bin_simulation, alpha = 1.0, normed = True, label = 'Simulation')
plt.hist(observation_data, color = 'k', bins = bin_observation, alpha = 0.5, normed = True,  label = 'Observation')
plt.legend()
plt.show()