.. _parameter:

Parameters
================

Universal Constants
---------------------

   .. py:function:: param.AU

      1 Astronomical Unit in cgs = 1.496e+13 

   .. py:function:: param.G

      Gravitational constant in cgs = 6.67e-8

   .. py:function:: param.M

      1 Solar Mass in cgs = 1.989e+33

   .. py:function:: param.M_H

      Mass of Hydrogen in cgs = 1.6733e-24

Initial Condition
--------------------

   Initial condition that can be changed. Here we adopted the case 7 from Visser (2009).

   .. py:function:: param.cs

      Sound Speed in cgs = 0.26*1e5

   .. py:function:: param.t_acc

      t_acc given in Visser (2009) = 2.5*1e5*31556926 [s]

   .. py:function:: param.omega

      Rotation rate = 1e-13

   .. py:function:: param.time

      The time to be calculated upon = [0.05*self.t_acc, 0.1*self.t_acc, 0.2*self.t_acc, 0.3*self.t_acc, 0.5*self.t_acc, 0.7*self.t_acc, self.t_acc]

   .. py:function:: param.r_env

      The outer edge of the simulation, where the TSC model breaks

   .. py:function:: param.max_tau

      The maximum tau value in TSC model, where to TSC perturbation breaks = 0.01

Variable for Calculation
-------------------------

   Formula for smapling points to solve on for simulation
   
   `z_g1_grid = [2/(i+1) for i in np.logspace(param.sample_range, 0, param.sample_pts)]`

   `z_s1_grid = np.logspace(0, param.sample_range, param.sample_pts)`

   .. py:function:: param.sample_pts

      Radial resolution inside and outside of the expansion wave front, separately = 1000

   .. py:function:: param.sample_range

      Sampling from 1e-x to 1ex = 4

   .. py:function:: param.theta_pts

      Angular resoliution calculated = 1000

Variable for Graphing
-----------------------

   .. py:function:: param.cgs_contour

      Contour value for cgs graph = 1e-18
   
   .. py:function:: param.cgs_vmin

      Minimum value for colormap for cgs graph = 1e-22

   .. py:function:: param.cgs_vmax

      Maximum value for colormap for cgs graph = 1e-
      
   .. py:function:: param.color_1d

      Color Value for different time graphed = ['blue', 'orange', 'green', 'red', 'purple', 'cyan', 'grey']

   .. py:function:: param.visser_contour

      Contour value for density of hydrogen graph = 1e6

   .. py:function:: param.visser_vmin

      Minimum value for colormap for density of hydrogen graph = 1e4

   .. py:function:: param.visser_vmax

      Maximum value for colormap for density of hydrogen graph = 1e9

.. automodule:: parameter
   :members:
   :undoc-members:
   :show-inheritance:
