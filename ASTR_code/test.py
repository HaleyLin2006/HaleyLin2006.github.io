# from parameter import param
# import matplotlib.pyplot as plt

# param = param()

# param.sample_pts = 10
# z_g1_grid = [2/(i+1) for i in np.logspace(param.sample_range, 0, param.sample_pts)]  # also x_s1
# z_s1_grid = np.logspace(0, param.sample_range, param.sample_pts)

# z_g1 = []
# z_s1 = []
# for i in range(param.sample_pts):
#     if i == 0:
#         z_s1_i = 0.99*z_s1_grid[i]
#         z_g1_i = 0.99*z_g1_grid[i]
#         z_s1.append(z_s1_i)
#         z_g1.append(z_g1_i)
#     else:
#         z_s1_i = 0.5*(z_s1_grid[i-1] + z_s1_grid[i])
#         z_g1_i = 0.5*(z_g1_grid[i-1] + z_g1_grid[i])
#         z_s1.append(z_sx1_i)
#         z_g1.append(z_g1_i)

# r_s1 = []
# r_g1 = []
# r_s1_grid = []
# r_g1_grid = []
# for t in param.time:
#     r_s1.append([1/i*param.cs*t for i in z_s1])
#     r_g1.append([1/i*param.cs*t for i in z_g1])
#     r_s1_grid.append([1/i*param.cs*t for i in z_s1_grid])
#     r_g1_grid.append([1/i*param.cs*t for i in z_g1_grid])

# r = []
# r_grid = []
# for i in range(len(param.time)):
#     r.append([j for j in r_s1[i]][::-1]+[j for j in r_g1[i]][::-1])
#     r_grid.append([j for j in r_s1_grid[i]][::-1]+[j for j in r_g1_grid[i]][::-1])

# plt.plot(range(param.sample_pts), [i for i in z_s1], '--bo')
# plt.plot(range(param.sample_pts), [i for i in z_s1_grid], '--bo', color='red')

# print(z_s1)
# print(z_s1_grid)

# # plt.xscale('log')
# plt.yscale('log')
# plt.show()

import numpy as np
import matplotlib.pyplot as plt
from parameter import param

param = param()

z_g1_grid = [2/(i+1) for i in np.logspace(param.sample_range, 0, param.sample_pts)]  # also x_s1
z_s1_grid = np.logspace(0, param.sample_range, param.sample_pts)

r_grid = []
for t in param.time:
    r_grid.append([1/i*param.cs*t for i in z_s1_grid[::-1]]+[1/i*param.cs*t for i in z_g1_grid[::-1]])

# for j in range(len(param.time)):
plt.plot(range(2*param.sample_pts), [i for i in r_grid[0]], '--bo', color='red')
plt.show()
