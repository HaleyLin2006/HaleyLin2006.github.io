.. _graph:
.. latex:samepage::

Graphing
===============

The graphing module consist of all graphing function for individual component and the entirety of simulation, 
combining from different paper. There are mainly two classes so far: `shu_graph` from Shu, 1977 is the envelope 
envolution without rotation, and `tsc_graph` from Terebey, Shu, Cassen (1984) is envelope evolution with small rotation.

You can change the initial condition given, such as time, initial mass, angular velocity..., in parameter_.

`shu_graph`
-------------
Initialize and graph Shu (1977), using Fluid model to model the gravitational collapse of isothermal cloud with no rotation.
   :Parameter(s): (`bool` or `None`) Optional. Determine whether to save the graph. Default is False.

   .. py:function:: shu_graph.dimensionless_a_v_1d

      Plot the dimensionless density and velocity profiles in 1D, against the radius (x=1/z).
      
      :param None: 

   .. py:function:: shu_graph.cgs_density_1d

      Plot the density profiles in CGS units in 1D for different time.

      :param None: 

   .. py:function:: shu_graph.nh_density_1d

      Plot the density profiles in density of hydrogen in 1D for different time.

      :param None: 

   .. py:function:: shu_graph.cgs_density_2d

      Plot density profiles in CGS units in 2D for different times using a logarithmic colormap with contour value (with `plt.imshow()`). 
      Note that it shows the entire range we calculated, therefore, the distance it shows is in log scale.

      :param None: 

   .. py:function:: shu_graph.nh_density_2d
      
      Plot density profiles in density of hydrogen in 2D for different times using 
      a logarithmic colormap with contour value (with `plt.contourf()`). Note that it only show the range in Visser (2009).

      :param None: 

`tsc_graph`
------------
   Initialize Terebey, Shu & Cassen (1984), using Fluid model to model the 
   gravitational collapse of isothermal cloud with small rotation.

   :Parameter(s): (`bool` or `None`)  Optional. Determine whether to save the graph. Default is False.

   .. py:function:: tsc_graph.dimensionless_a_1d

      Plot the dimensionless density profiles in 1D against radius (1/z).

      :param None: 
   
   .. py:function:: tsc_graph.dimensionless_quad_v_w_1d

      Plot dimensionless quadrupolar velocity in theta and phi direction profiles in 1D against radius (x=1/z).

      :param None: 

   .. py:function:: tsc_graph.dimensionless_mono_v_1d

      Plot dimensionless monopolar velocity profiles in 1D against radius (x=1/z).

      :param None: 

   .. py:function:: tsc_graph.cgs_density_1d

      Plot density profiles in CGS units in 1D for different times.

      :param None: 
   
   .. py:function:: tsc_graph.nh_density_1d

      Plot density profiles in density of hydrogen in 1D for different times.

      :param None: 

   .. py:function:: tsc_graph.cgs_density_2d

      Plot density profiles in CGS units in 2D for different times using a logarithmic colormap 
      with contour value (with `plt.imshow()`). Note that it shows the entire range we calculated, therefore, the distance it shows is in log scale.

      :param None: 
   
   .. py:function:: tsc_graph.nh_density_2d

      Plot density profiles in density of hydrogen units in 2D for different times using 
      contour plots with a logarithmic colormap (with `plt.countourf()`). Note that it only show the range in Visser (2009).

      :param None: 

   .. py:function:: tsc_graph.int_mass

      Integrate the mass of envelope from inner edge to outer edge, to keep track of the mass of the disk, and make sure it fits the initial condition.

      :param None: 
      :return: The mass of disk in cgs.



.. automodule:: graphing
   :members:
   :undoc-members:
   :show-inheritance:
