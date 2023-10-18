import numpy as np

class param():
    def __init__(self):
        #   universal constant in cgs
        self.AU     = 1.496e+13                            # 1 Astronomical Unit            [cm]
        self.G      = 6.67e-8                              # Gravitational Constant       
        self.M      = 1.989e+33                            # 1 solar mass                   [g]
        self.M_H    = 1.6733e-24                           # Mass of hydrogen               [g]

        #   variable for model
        self.cs             = 0.26*1e5                                                              # sound speed in Visser (2009)   [cm/s]
        self.t_acc          = 2.5*1e5 * 31556926                                                    # t_acc in Visser (2009)         [s]
        self.omega          = 1e-13                                                                 # rotation rate in Visser (2009) [s-1]
        self.r_env          = self.G*self.M/2/self.cs**2 / self.AU                                  # outer radius 
        self.time           = [0.05*self.t_acc, 0.1*self.t_acc, 0.2*self.t_acc, 0.3*self.t_acc, 
                               0.5*self.t_acc, 0.7*self.t_acc, self.t_acc]                          # time in Visser (2009) [s]
        # self.time           = [0.1*self.t_acc, 0.5*self.t_acc, self.t_acc]
        self.max_tau        = 0.01

        #   varaible for calculation
        self.sample_pts     = 1000
        self.sample_range   = 4
        self.theta_pts      = 1000

        #   variable for graph
        self.cgs_contour    = 1e-18
        self.cgs_vmin       = 1e-22
        self.cgs_vmax       = 1e-14
        self.color_1d       = ['blue', 'orange', 'green', 'red', 'purple', 'cyan', 'grey']

        self.visser_contour = 1e6
        self.visser_vmin    = 1e4
        self.visser_vmax    = 1e9

        self.tick           = [0.0, 0.5, 1.0, 1.5]