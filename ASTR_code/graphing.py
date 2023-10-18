import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
import numpy as np

from envelope.shu77 import shu
from envelope.tsc84 import tsc
from data_handler import data
from parameter import param

# TODO: save data using npz
#       -> create folder for each different param + a csv record param
# TODO: sphinx

# TODO: integrate and check total mass

param = param()
data = data()

z_g1_grid = [2/(i+1) for i in np.logspace(param.sample_range, 0, param.sample_pts)]  # also x_s1
z_s1_grid = np.logspace(0, param.sample_range, param.sample_pts)

r_grid = []
for t in param.time:
    r_grid.append([1/i*param.cs*t for i in z_s1_grid[::-1]]+[1/i*param.cs*t for i in z_g1_grid[::-1]])

z_g1 = []
z_s1 = []
for i in range(param.sample_pts):
    if i == 0:
        z_s1_i = 0.5*z_s1_grid[i]
        z_g1_i = 0.5*z_g1_grid[i]
        z_s1.append(z_s1_i)
        z_g1.append(z_g1_i)
    else:
        z_s1_i = 0.5*(z_s1_grid[i-1] + z_s1_grid[i])
        z_g1_i = 0.5*(z_g1_grid[i-1] + z_g1_grid[i])
        z_s1.append(z_s1_i)
        z_g1.append(z_g1_i)

#   Envelope
class shu_graph():
    def __init__(self, save=False):
        """
        Initialize Shu (1977), using Fluid model to model the 
        gravitational collapse of isothermal cloud with no rotation.

        Args:
            save (bool): Whether to save the graph. Default is False.

        (Output):
            self.radius     [cm]/AU
            self.nh_density [nh] -> *param.M_H to go back to cgs
            self.v_r        [cgs]
            
        """
        data_param = [param.cs, param.t_acc, param.omega, param.time]
        if data.isSaved(param=data_param, model="Shu") == True:
            self.radius_str, self.nh_density_str, self.v_r_str, self.v_theta_str, self.v_phi_str = data.load(param=data_param, model="Shu")
            
            # print(len(self.radius), len(self.radius[0]))
            self.radius     = [list(map(float, line.split('\t'))) for line in self.radius_str.split('|')]
            self.nh_density = [list(map(float, line.split('\t'))) for line in self.nh_density_str.split('|')]
            self.v_r        = [list(map(float, line.split('\t'))) for line in self.v_r_str.split('|')]

            print(len(self.radius), len(self.radius[0]))
            print(len(self.nh_density), len(self.nh_density[0]))
            print("Loading Successfully!")

        if data.isSaved(param=data_param, model="Shu") == False:
            # Create a Shu instance with specified parameters
            self.shu = shu(z_g1=z_g1, z_s1=z_s1, cs=param.cs, time=param.time)

            # Initialize self-similar variables
            self.z_g1       = self.shu.get_z_g1()
            self.z_s1       = self.shu.get_z_s1()
            self.a_shu_s1   = self.shu.get_a_s1()
            self.v_shu_s1   = self.shu.get_v_s1()
            self.a_shu_g1   = self.shu.get_a_g1()
            self.v_shu_g1   = self.shu.get_v_g1()
            self.dadz       = self.shu.get_dadz()
            self.dvdz       = self.shu.get_dvdz()

            # Initialize dimension variables
            self.radius_g1      = self.shu.get_r_g1()
            self.radius_s1      = self.shu.get_r_s1()
            self.cgs_density_g1 = self.shu.get_cgs_rho_g1()
            self.cgs_density_s1 = self.shu.get_cgs_rho_s1()
            self.nh_density_g1  = self.shu.get_nh_rho_g1()
            self.nh_density_s1  = self.shu.get_nh_rho_s1()

            self.radius     = []
            self.nh_density = []
            self.v_r        = []
            for i in range(len(param.time)):
                self.radius.append([j/param.AU for j in self.radius_s1[i]][::-1]+[j/param.AU for j in self.radius_g1[i]][::-1])
                self.nh_density.append(self.nh_density_s1[i][::-1]+self.nh_density_g1[i][::-1])
                self.v_r.append([param.cs*j for j in self.v_shu_s1][::-1]+[param.cs*j for j in self.v_shu_g1][::-1])

            if save == True and data.isSaved(param=data_param, model="Shu") == False:

                self.radius_str     = '|'.join('\t'.join(map(str, sublist)) for sublist in self.radius)
                self.nh_density_str = '|'.join('\t'.join(map(str, sublist)) for sublist in self.nh_density)
                self.v_r_str        = '|'.join('\t'.join(map(str, sublist)) for sublist in self.v_r)

                data.save(param=data_param, model="Shu", radius=self.radius_str, density=self.nh_density_str, v_r=self.v_r_str)
    
                self.radius     = [list(map(float, line.split('\t'))) for line in self.radius_str.split('|')]
                self.nh_density = [list(map(float, line.split('\t'))) for line in self.nh_density_str.split('|')]
                self.v_r        = [list(map(float, line.split('\t'))) for line in self.v_r_str.split('|')]

                print(len(self.radius), len(self.radius[0]))
                print(len(self.nh_density), len(self.nh_density[0]))

    def dimensionless_a_v_1d(self):
        """
        Plot the dimensionless density and velocity profiles in 1D,
        against the radius (1/z).
        """
        self.shu = shu(z_g1=z_g1, z_s1=z_s1, cs=param.cs, time=param.time)

        self.z_g1       = self.shu.get_z_g1()
        self.z_s1       = self.shu.get_z_s1()
        self.a_shu_s1   = self.shu.get_a_s1()
        self.v_shu_s1   = self.shu.get_v_s1()
        self.a_shu_g1   = self.shu.get_a_g1()
        self.v_shu_g1   = self.shu.get_v_g1()

        plt.plot([1/i for i in self.z_g1], self.a_shu_g1)
        plt.plot([1/i for i in self.z_s1], self.a_shu_s1)
        plt.plot([1/i for i in self.z_g1], [-i for i in self.v_shu_g1])
        plt.plot([1/i for i in self.z_s1], [-i for i in self.v_shu_s1])

    def cgs_density_1d(self):
        """
        Plot the density profiles in CGS units in 1D for different time.
        """
        for i in range(len(param.time)):
            plt.plot(self.radius[i], [j*param.M_H for j in self.nh_density[i]], color=param.color_1d[i], label=f't={param.time[i]/1e12}e12')

    def nh_density_1d(self):
        """
        Plot the density profiles in density of hydrogen in 1D for different time.
        """
        for i in range(len(param.time)):
            plt.plot(self.radius[i], self.nh_density[i], color=param.color_1d[i], label=f't={param.time[i]/1e12}e12')

    def cgs_density_2d(self):
        """
        Plot density profiles in CGS units in 2D for different times using 
        a logarithmic colormap with contour value. Note that it shows the entire range we calculated

        LOG SCALE!
        """
        fig, axes = plt.subplots(2, 4)
        axes_flat = axes.flatten()

        for i, ax in enumerate(axes_flat):
            if i < len(param.time):
                print(f"progress: {i+1}/{len(param.time)}, currently at {param.time[i]/param.t_acc} t_acc")
                
                radii = np.linspace(0, 2*param.sample_pts, 2*param.sample_pts)  # Distances from the origin
                molecule_counts = [j*param.M_H for j in self.nh_density[i]]

                # Create a grid
                x = np.linspace(0, 2*param.sample_pts, 2*param.sample_pts)
                y = np.linspace(0, 2*param.sample_pts, 2*param.sample_pts)
                xx, yy = np.meshgrid(x, y)
                r = np.sqrt(xx**2 + yy**2)

                # Interpolate data onto the grid
                interpolated_data = np.interp(r.ravel(), radii, molecule_counts, left=0, right=0)
                interpolated_data = interpolated_data.reshape(xx.shape)

                # Logarithmic colormap
                norm = LogNorm(vmin=param.cgs_vmin, vmax=param.cgs_vmax)

                # Plot the radially symmetrical image
                im = ax.imshow(interpolated_data, cmap='gist_rainbow_r', origin='lower', 
                               extent=[0, 2*param.sample_pts, 0, 2*param.sample_pts], norm=norm)

                # Add contour line at the value 10^(-18)
                contour_value = param.cgs_contour 
                ax.contour(xx, yy, interpolated_data, levels=[contour_value], colors='red', linewidths=2)

                ax.set_title(f't = {param.time[i]/param.t_acc}')

            else:
                plt.axis("off")
                cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
                cbar.set_label('Density')

        fig.suptitle('Shu (1977) 2D Graph (cgs)')

    def nh_density_2d(self):
        """
        Plot density profiles in density of hydrogen in 2D for different times using 
        a logarithmic colormap with contour value. Note that it only show the range in Visser.
        """

        fig, axes = plt.subplots(2, 2)
        axes_flat = axes.flatten()

        for i, ax in enumerate(axes_flat):
            if i < len(param.time):
                print(f"progress: {i+1}/{len(param.time)}, currently at {param.time[i]/param.t_acc} t_acc")
                radii = [j/1e3 for j in self.radius[i]]  # Distances from the origin
                molecule_counts = self.nh_density[i]
                
                # Create a grid
                x = [j/1e3 for j in self.radius[i]] 
                y = [j/1e3 for j in self.radius[i]] 

                xx, yy = np.meshgrid(x, y)
                r = np.sqrt(xx**2 + yy**2)

                # Interpolate data onto the grid
                interpolated_data = np.interp(r.ravel(), radii, molecule_counts, left=0, right=0)
                interpolated_data = interpolated_data.reshape(xx.shape)

                # Logarithmic colormap
                norm = LogNorm(vmin=param.visser_vmin, vmax=param.visser_vmax)

                # Plot the radially symmetrical image
                im = ax.contourf(x, y, interpolated_data, cmap='gist_rainbow_r', levels=np.logspace(4,9,100), norm=norm, aspect='auto')
                
                ax.set_xlim(0, 1.5)
                ax.set_ylim(0, 1.5)

                r_env_edge = patches.Circle((0,0), param.r_env/1e3, fill=False, edgecolor='black', linewidth=2)
                ax.add_patch(r_env_edge)

                # Add contour line at the value
                ax.contour(xx, yy, interpolated_data, levels=[param.visser_contour],         
                           colors='red', linewidths=2)

                ax.set_title(f't = {param.time[i]/param.t_acc}')

            else:
                plt.axis("off")
                cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
                cbar.set_ticks([10**4, 10**5, 10**6, 10**7, 10**8, 10**9])
                cbar.set_label('Density')

        fig.suptitle('Shu (1977) 2D Graph (n_h)')

class tsc_graph():
    def __init__(self, save=False):
        """
        Initialize Terebey, Shu & Cassen (1984), using Fluid model to model the 
        gravitational collapse of isothermal cloud with small rotation.

        Args:
            save (bool): Whether to save the graph. Default is False.

        (Output):
            self.radius     [cm]/AU
            self.nh_density [nh] -> *param.M_H to go back to cgs
            self.v_r        [cgs]
        """
        data_param = [param.cs, param.t_acc, param.omega, param.time]
        if data.isSaved(param=data_param, model="TSC") == True:
            self.radius_str, self.nh_density_str, self.v_r_str, self.v_theta_str, self.v_phi_str = data.load(param=data_param, model="TSC")

            # Convert the string back to a nested list
            self.radius         = [list(map(float, line.split('\t'))) for line in self.radius_str.split('|')]
            self.nh_density     = [[[float(num) for num in sublist.split('\t')] for sublist in part.split(';')] for part in self.nh_density_str.split('|')]
            self.cgs_v_r        = [[[float(num) for num in sublist.split('\t')] for sublist in part.split(';')] for part in self.v_r_str.split('|')]
            self.cgs_v_theta    = [[[float(num) for num in sublist.split('\t')] for sublist in part.split(';')] for part in self.v_theta_str.split('|')]
            self.cgs_v_phi      = [[[float(num) for num in sublist.split('\t')] for sublist in part.split(';')] for part in self.v_phi_str.split('|')]
            
            print("radius",     len(self.radius),       len(self.radius[0]))
            print("density",    len(self.nh_density),   len(self.nh_density[0]))
            print("v_r",        len(self.cgs_v_r),      len(self.cgs_v_r[0]))
            print("v_theta",    len(self.cgs_v_theta),  len(self.cgs_v_theta[0]))
            print("v_phi",      len(self.cgs_v_phi),    len(self.cgs_v_phi[0]))
            
            print("Loading Successfully!")
        
        if data.isSaved(param=data_param, model="TSC") == False:
            # Create Shu and TSC instances with specified parameters
            self.shu = shu(z_g1=z_g1, z_s1=z_s1, cs=param.cs, time=param.time)
            self.tsc = tsc(z_g1=z_g1, z_s1=z_s1, a_shu_s1=self.shu.get_a_s1(), a_shu_g1=self.shu.get_a_g1(), 
                            v_shu_s1=self.shu.get_v_s1(), v_shu_g1=self.shu.get_v_g1(),
                            store_dadz=self.shu.get_dadz(), store_dvdz=self.shu.get_dvdz(),
                            cs=param.cs, time=param.time)

            # Initialize dimension variables
            self.radius_s1      = self.tsc.get_r_s1()
            self.radius_g1      = self.tsc.get_r_g1()

            self.cgs_density_s1 = self.tsc.get_cgs_density_s1()
            self.cgs_density_g1 = self.tsc.get_cgs_density_g1()
            self.nh_density_s1  = self.tsc.get_nh_density_s1()
            self.nh_density_g1  = self.tsc.get_nh_density_g1()

            self.cgs_v_r_s1     = self.tsc.get_cgs_v_r_s1()
            self.cgs_v_r_g1     = self.tsc.get_cgs_v_r_g1()
            self.cgs_v_theta_s1 = self.tsc.get_cgs_v_theta_s1()
            self.cgs_v_theta_g1 = self.tsc.get_cgs_v_theta_g1()
            self.cgs_v_phi      = self.tsc.get_cgs_v_phi()

            self.inner_edge     = [(param.omega*t)**2* param.cs*t/param.AU for t in param.time]
            self.outer_edge     = [1/(param.omega*t) * param.cs*t/param.AU for t in param.time]
            print(param.omega*param.M*param.G/0.975/param.cs**3)
            print([param.omega*t for t in param.time])
            print(self.outer_edge)

            self.radius         = []
            self.nh_density     = []
            self.cgs_v_r        = []
            self.cgs_v_theta    = []
            self.cgs_v_phi      = self.cgs_v_phi

            for i in range(len(param.time)):
                print(f"reshape process, {i}/{len(param.time)}" )
                self.radius.append([j/param.AU for j in self.radius_s1[i]][::-1]+[j/param.AU for j in self.radius_g1[i]][::-1])

                self.nh_density.append([])
                self.cgs_v_r.append([])
                self.cgs_v_theta.append([])
                for j in range(param.theta_pts):
                    self.nh_density[-1].append(self.nh_density_s1[i][j][::-1]+self.nh_density_g1[i][j][::-1])
                    self.cgs_v_r[-1].append(self.cgs_v_r_s1[i][j][::-1]+self.cgs_v_r_g1[i][j][::-1])
                    self.cgs_v_theta[-1].append(self.cgs_v_theta_s1[i][j][::-1]+self.cgs_v_theta_g1[i][j][::-1])
                
            print(len(r_grid[0]), len(self.radius[0]))

            print("radius",     len(self.radius),       len(self.radius[0]))
            print("density",    len(self.nh_density),   len(self.nh_density[0]))
            print("v_r",        len(self.cgs_v_r),      len(self.cgs_v_r[0]))
            print("v_theta",    len(self.cgs_v_theta),  len(self.cgs_v_theta[0]))
            print("v_phi",      len(self.cgs_v_phi),    len(self.cgs_v_phi[0]))
            
            if save == True and data.isSaved(param=data_param, model="TSC") == False:
                
                self.radius_str         = '|'.join('\t'.join(map(str, sublist)) for sublist in self.radius)
                self.nh_density_str     = '|'.join(';'.join('\t'.join(map(str, inner_list)) for inner_list in middle_list) for middle_list in self.nh_density)
                self.cgs_v_r_str        = '|'.join(';'.join('\t'.join(map(str, inner_list)) for inner_list in middle_list) for middle_list in self.cgs_v_r)
                self.cgs_v_theta_str    = '|'.join(';'.join('\t'.join(map(str, inner_list)) for inner_list in middle_list) for middle_list in self.cgs_v_theta)
                self.cgs_v_phi_str      = '|'.join(';'.join('\t'.join(map(str, inner_list)) for inner_list in middle_list) for middle_list in self.cgs_v_phi)

                data.save(param=data_param, model="TSC", radius=self.radius_str, density=self.nh_density_str, 
                          v_r=self.cgs_v_r_str, v_theta=self.cgs_v_theta_str, v_phi=self.cgs_v_phi_str)
                print("Saving Successfully!")

                self.radius         = [list(map(float, line.split('\t'))) for line in self.radius_str.split('|')]
                self.nh_density     = [[[float(num) for num in sublist.split('\t')] for sublist in part.split(';')] for part in self.nh_density_str.split('|')]
                self.cgs_v_r        = [[[float(num) for num in sublist.split('\t')] for sublist in part.split(';')] for part in self.cgs_v_r_str.split('|')]
                self.cgs_v_theta    = [[[float(num) for num in sublist.split('\t')] for sublist in part.split(';')] for part in self.cgs_v_theta_str.split('|')]
                self.cgs_v_phi      = [[[float(num) for num in sublist.split('\t')] for sublist in part.split(';')] for part in self.cgs_v_phi_str.split('|')]

            print("radius",     len(self.radius),       len(self.radius[0]))
            print("density",    len(self.nh_density),   len(self.nh_density[0]))
            print("v_r",        len(self.cgs_v_r),      len(self.cgs_v_r[0]))
            print("v_theta",    len(self.cgs_v_theta),  len(self.cgs_v_theta[0]))
            print("v_phi",      len(self.cgs_v_phi),    len(self.cgs_v_phi[0]))
       
    def dimensionless_a_1d(self):
        """
        Plot the dimensionless density profiles in 1D against radius (1/z).
        """
        # Plot Shu acceleration profiles
        plt.plot([1/i for i in z_g1], [i for i in self.shu.get_a_g1()], color='orange', label='a shu', linewidth=2)
        plt.plot([1/i for i in z_s1], [i for i in self.shu.get_a_s1()], color='orange', linewidth=2)
        
        # Plot quadrupolar acceleration profiles
        plt.plot([1/i for i in z_g1], [i for i in self.quad_a_g1], label=r"$-\alpha_Q$", color='red', linewidth=2)
        plt.plot([1/i for i in z_s1], [abs(i) for i in self.quad_a_s1], color='red', linewidth=2)                 
        
        # Plot monopolar acceleration profiles
        plt.plot([1/i for i in z_g1], [i for i in self.mono_a_g1], label=r"$a_m$", color='blue', linewidth=2)
        plt.plot([1/i for i in z_s1], [i for i in self.mono_a_s1], color='blue', linewidth=2)
    
    def dimensionless_quad_v_w_1d(self):
        """
        Plot dimensionless quadrupolar velocity in theta and phi 
        direction profiles in 1D against radius (1/z).
        """
        # Plot quadrupolar velocity profiles
        plt.plot([1/i for i in z_g1], [i for i in self.quad_v_g1], label=r"$v_Q$", color='blue', linewidth=2)
        plt.plot([1/i for i in z_s1], [i for i in self.quad_v_s1], color='blue', linewidth=2)       
        
        # Plot quadrupolar angular velocity profiles
        plt.plot([1/i for i in z_g1], [i for i in self.quad_w_g1], label=r"$-w_Q$", color='green', linewidth=2)
        plt.plot([1/i for i in z_s1], [-i for i in self.quad_w_s1], color='green', linewidth=2)

    def dimensionless_mono_v_1d(self):
        """
        Plot dimensionless monopolar velocity profiles in 1D against radius (1/z).
        """
        # Plot Shu and monopolar velocity profiles
        plt.plot([1/i for i in z_s1], [-i for i in self.shu.get_v_s1()], color='orange', label='-v shu', linewidth=2)
        plt.plot([1/i for i in z_s1], [i for i in self.mono_v_s1], label=r"$v_m$", linewidth=2)

    def cgs_density_1d(self):
        """
        Plot density profiles in CGS units in 1D for different times.
        """
        for i in range(len(param.time)):
            plt.plot(self.radius[i], [j*param.M_H for j in self.nh_density[i][0]], color=param.color_1d[i], label=f't={param.time[i]/1e12}e12')
    
    def nh_density_1d(self):
        """
        Plot density profiles in density of hydrogen in 1D for different times.
        """
        for i in range(len(param.time)):
            plt.plot(self.radius[i], self.nh_density[i][0], color=param.color_1d[i], label=f't={param.time[i]/1e12}e12')

    def cgs_density_2d(self):
        """
        Plot density profiles in CGS units in 2D for different times using 
        a logarithmic colormap with contour value.
        """
        fig, axes = plt.subplots(2, 2)
        axes_flat = axes.flatten()

        theta_values = np.linspace(0, 2 * np.pi, param.theta_pts)

        i = 0
        for i, ax in enumerate(axes_flat):
            if i < len(param.time):
                print(f"progress: {i+1}/{len(param.time)}, currently at {param.time[i]/param.t_acc} t_acc")
                data_values = [] 
                for theta in range(len(theta_values)):
                    data_values.append([j*param.M_H for j in self.nh_density[i][theta]])

                # Example data in (r, theta) format
                r_values = np.linspace(0, 2*param.sample_pts, 2*param.sample_pts)
                # data_values = [i*param.M_H for i in self.nh_density]  # List of lists, each sublist represents data at a theta value

                # Convert (r, theta) to (x, y) coordinates
                x_values = []
                y_values = []
                for j, theta in enumerate(theta_values):
                    x_values.extend(r_values * np.cos(theta))
                    y_values.extend(r_values * np.sin(theta))

                # Create a flattened (x, y) array and corresponding flattened data array
                x_flat = np.array(x_values)
                y_flat = np.array(y_values)
                data_flat = np.array(data_values).flatten()

                # Create a grid in (x, y) coordinates
                x_grid = np.linspace(0, max(x_flat), 2*param.sample_pts)
                y_grid = np.linspace(0, max(y_flat), 2*param.sample_pts)
                xx, yy = np.meshgrid(x_grid, y_grid)

                # Interpolate data onto the (x, y) grid
                data_interp = griddata((x_flat, y_flat), data_flat, (xx, yy), method='linear')

                # Plot the interpolated data using imshow()
                norm = LogNorm(vmin=param.cgs_vmin, vmax=param.cgs_vmax)
                im = ax.imshow(data_interp, extent=[0, max(x_flat), 0, max(y_flat)],
                        cmap='gist_rainbow_r', origin='lower', aspect='auto', norm=norm)
                
                contour_value = param.cgs_contour
                ax.contour(xx, yy, data_interp, levels=[contour_value], colors='red', linewidths=2)

                ax.set_title(f't = {param.time[i]/param.t_acc}')

            else:
                plt.axis('off')
                cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
                cbar.set_label('Density')

        fig.suptitle('TSC (1984) 2D Graph (cgs)')

    def nh_density_2d(self):
        """
        Plot density profiles in density of hydrogen units in 2D for different times using 
        contour plots with a logarithmic colormap.
        """
        fig, axes = plt.subplots(2, 4)
        axes_flat = axes.flatten()

        theta_values = np.linspace(0, 2*np.pi, param.theta_pts)

        i = 0
        for i, ax in enumerate(axes_flat):
            if i < len(param.time):
                print(f"progress: {i+1}/{len(param.time)}, currently at {param.time[i]/param.t_acc} t_acc")
                
                molecule_counts = [] 
                for theta in range(len(theta_values)):
                    molecule_counts.append(self.nh_density[i][theta])

                # Example data in (r, theta) format
                r_values =  np.array([j/1e3 for j in self.radius[i]])
                
                # Convert (r, theta) to (x, y) coordinates
                x_values = []
                y_values = []
                for j, theta in enumerate(theta_values):
                    x_values.extend(r_values * np.cos(theta))
                    y_values.extend(r_values * np.sin(theta))

                # Create a flattened (x, y) array and corresponding flattened data array
                x_flat = np.array(x_values)
                y_flat = np.array(y_values)
                data_flat = np.array(molecule_counts).flatten()

                # Create a grid in (x, y) coordinates
                x_grid = [j/1e3 for j in self.radius[i]]
                y_grid = [j/1e3 for j in self.radius[i]]
                xx, yy = np.meshgrid(x_grid, y_grid)

                # Interpolate data onto the (x, y) grid
                data_interp = griddata((x_flat, y_flat), data_flat, (xx, yy), method='linear')

                # Create a contour plot with a logarithmic colormap
                norm = LogNorm(vmin=param.visser_vmin, vmax=param.visser_vmax)
                im = ax.contourf(xx, yy, data_interp, cmap='gist_rainbow_r', levels=np.logspace(4,9,100), norm=norm)

                self.inner_edge     = [(param.omega*t)**2* param.cs*t/param.AU for t in param.time]
                self.outer_edge     = [1/(param.omega*t) * param.cs*t/param.AU for t in param.time]

                ax.set_xlim(0, 10)
                ax.set_ylim(0, 10)

                # Add contour line at the value
                ax.contour(xx, yy, data_interp, levels=[param.visser_contour],           
                           colors='red', linewidths=2)
                
                r_env_edge = patches.Circle((0,0), param.r_env/1e3, fill=False, edgecolor='white', linewidth=2)
                ax.add_patch(r_env_edge)
                circ_inner = patches.Circle((0,0), self.inner_edge[i]/1e3, fill=False, edgecolor='black', linewidth=2)
                ax.add_patch(circ_inner)
                circ_outer = patches.Circle((0,0), self.outer_edge[i], fill=False, edgecolor='black', linewidth=2)
                ax.add_patch(circ_outer)

                ax.set_title(f't = {param.time[i]/param.t_acc}')
                
            else:
                plt.axis('off')
                cbar = plt.colorbar(im, ax=axes.ravel().tolist(), orientation='vertical')
                cbar.set_ticks([10**4, 10**5, 10**6, 10**7, 10**8, 10**9])
                cbar.set_label('Density')
    
        fig.suptitle('TSC (1984) 2D Graph (n_h)')
    
    def int_mass(self):
        # 2*pi S(rho * r^2)dr
        for i in range(len(self.time)):
            for j in range(len(self.radius)):
                self.radius
                r_grid
        
        return     

########## Shu (1977) Graph ##########
# s = shu_graph(save=True)
# s.dimensionless_a_v_1d()
# s.cgs_density_1d()
# s.nh_density_1d()
# s.cgs_density_2d()
# s.nh_density_2d()

########## TSC (1984) Graph ##########
t = tsc_graph(save=False)
# t.dimensionless_a_1d()
# t.dimensionless_quad_v_w_1d()
# t.dimensionless_mono_v_1d()
# t.cgs_density_1d()
# t.nh_density_1d()
# t.cgs_density_2d()
# t.int_mass()
t.nh_density_2d()

plt.xscale('log')
plt.yscale('log')
plt.show()