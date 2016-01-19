import json

rsets = dict()
rsets['half length'] = 150. # [m] half length of model so that xmin/ymin = -half length and xmax/ymax = half length
rsets['density'] = 0.0001#0.0015 # features per square meter
rsets['shear angle'] = 10. # [deg]
rsets['rmin'] = 1. # [m] minimum dome radius
rsets['rmax'] = 2. # [m] maximum dome radius
rsets['rsigma'] = 1.5 # [m] standard deviation of radii distribution
rsets['orientation sigma'] = 10. # [deg] standard deviation of orientation distribuition
rsets['orientation mean'] = 0. # [deg] mean orientation of features w.r.t. x-axis
rsets['walls'] = 10 # number of double walls east and north
rsets['opposing walls distance'] = 7. # [m] distance between two opposing walls
rsets['height scaling factor'] = 0.29 # 0-1, where 1 means no scaling
rsets['dome'] = False #True # True for domes, False for ridges
rsets['ridges aspect ratio'] = 20.
rsets['ridges bending angle'] = 30.

with open('rhino_settings.json', 'w') as f:
    f.write(json.dumps(rsets, indent=2, sort_keys=True))