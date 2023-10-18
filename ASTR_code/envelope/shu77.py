from scipy.integrate import ode
import numpy as np
import sys

sys.path.insert(1, '/Users/haleylin/Desktop/astro')
from parameter import param

param = param()

class shu():
    # solve the EOM with z=1/x parameter provided by Shu (1977, S77)
    #   z_g1  (-)   : a list of radius the EOM will be solve on for z>1
    #   z_s1  (-)   : a list of radius the EOM will be solve on for z<1
    #   cs ([cm/s]) : sound speed given by S77 to bring dimensionless solution back with dimension
    #   time ([s])  : time is provided to bring dimensionless solution back with dimension
    #   num_data_pts: number of data points used
    def __init__(self, z_g1, z_s1, cs=0.2*1e5, time=[1e12, 2e12, 4e12, 8e12]):
        self.z_g1 = z_g1
        self.z_s1 = z_s1
        self.cs = cs
        self.t = time

        self.store_dadz = []
        self.store_dvdz = []
        def shu_ode(z, param):
            a = param[0]
            v = param[1]

            dadz = -a/((1-z*v)**2 - z**2) * ((1/z-v) * (a - 2*z*(1/z-v)))
            dvdz = -1/((1-z*v)**2 - z**2) * ((1/z-v) * (a*(1/z-v) - 2*z))

            if dadz<0:
                dadz = 0

            #global store_dadz
            self.store_dadz.append(dadz)
            #global store_dvdz
            self.store_dvdz.append(dvdz)

            # progress recorder
            print("shu, z:", z)

            return [dadz, dvdz]

        # initial condition for Shu, z -> 0 (x -> infty)
        A = 2.0001
        param0_g1 = [A*(self.z_g1[0]**2), -(A-2)*self.z_g1[0]]

        # Solution for Shu, z<1, (x>1)
        sol_shu_g1 = np.zeros((len(self.z_g1), len(param0_g1)))
        sol_shu_g1[0, :] = param0_g1
        f_shu_g1 = ode(shu_ode)
        f_shu_g1.set_integrator("vode", method = "bdf", atol=1e-20, nsteps=1000)
        f_shu_g1.set_initial_value(param0_g1, self.z_g1[0])

        for i in range(1, len(self.z_g1)):
            sol_shu_g1[i, :] = f_shu_g1.integrate(self.z_g1[i])
            f_shu_g1.set_initial_value(sol_shu_g1[i, :], self.z_g1[i])

        self.a_shu_g1 = [i for i in sol_shu_g1[:, 0]]
        self.v_shu_g1 = [i for i in sol_shu_g1[:, 1]]

        self.store_dadz = []
        self.store_dvdz = []

        # initial condition for Shu, z>1, (x<1)
        param0_s1 = [self.a_shu_g1[-1], self.v_shu_g1[-1]]

        # Solution for Shu, z>1, (x<1)
        sol_shu_s1 = np.zeros((len(self.z_s1), len(param0_s1)))
        sol_shu_s1[0, :] = param0_s1
        f_shu_s1 = ode(shu_ode)
        f_shu_s1.set_integrator("vode", method="bdf", atol=1e-20, nsteps=1000)
        f_shu_s1.set_initial_value(param0_s1, self.z_s1[0])

        for i in range(1, len(self.z_s1)):
            sol_shu_s1[i, :] = f_shu_s1.integrate(self.z_s1[i])
            f_shu_s1.set_initial_value(sol_shu_s1[i, :], self.z_s1[i])

        self.a_shu_s1 = [i for i in sol_shu_s1[:, 0]]
        self.v_shu_s1 = [i for i in sol_shu_s1[:, 1]]

    # return dadz x<1
    def get_dadz(self):
        return self.store_dadz
    
    # return dvdz x<1
    def get_dvdz(self):
        return self.store_dvdz

    ###### Radius #####
    # return dimensionless radius x>1
    def get_z_g1(self):
        return self.z_g1
    # return radius with dimension [cm] for x>1
    def get_r_g1(self):
        r = []
        for t in self.t:
            r.append([1/i*self.cs*t for i in self.z_g1])
        
        return r
    
    # return dimensionless radius x<1
    def get_z_s1(self):
        return self.z_s1
    # return radius with dimension [cm] for x<1
    def get_r_s1(self):
        r = []
        for t in self.t:
            r.append([1/i*self.cs*t for i in self.z_s1])
        
        return r

    ###### Density #####
    # return dimensionless parametric solution: density x>1
    def get_a_g1(self):
        return self.a_shu_g1
    
    # return density with dimension [g/cm^3] x>1
    def get_cgs_rho_g1(self):
        rho = []
        for t in self.t:
            rho.append([(i/(4*np.pi*param.G*t**2)) for i in self.a_shu_g1])
        
        return rho
    
    # return density with dimension of density of hydrogen x>1
    def get_nh_rho_g1(self):
        rho = []
        for t in self.t:
            rho.append([(i/(4*np.pi*param.G*t**2))/param.M_H for i in self.a_shu_g1])
        
        return rho

    # return dimensionless parametric solution: density x<1
    def get_a_s1(self):
        return self.a_shu_s1
    
    # return density with dimension [g/cm^3] x<1
    def get_cgs_rho_s1(self):
        rho = []
        for t in self.t:
            rho.append([(i/(4*np.pi*param.G*t**2)) for i in self.a_shu_s1])
        
        return rho
    
    def get_nh_rho_s1(self):
        rho = []
        for t in self.t:
            rho.append([(i/(4*np.pi*param.G*t**2))/param.M_H for i in self.a_shu_s1])
        
        return rho
    
    ###### Radial Velocity #####
    # return dimensionless pramatric solution: velocity x>1
    def get_v_g1(self):
        return self.v_shu_g1
    # return velocity with dimension [cm/s] x>1
    def get_u_g1(self):
        u = []
        for t in self.t:
            u.append([-(self.cs*i) for i in self.v_shu_g1])

        return u

    # return dimensionless pramatric solution: velocity x<1
    def get_v_s1(self):
        return self.v_shu_s1
    # return velocity with dimension [cm/s] x<1
    def get_u_s1(self):
        u = []
        for t in self.t:
            u.append([-(self.cs*i) for i in self.v_shu_s1])

        return u