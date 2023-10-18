from scipy.integrate import odeint, ode
from scipy.special import legendre
from numpy.polynomial import Polynomial
import numpy as np
import sys

sys.path.insert(1, '/Users/haleylin/Desktop/astro')
from parameter import param

param = param()

class tsc():
    def __init__(self, z_g1, z_s1, a_shu_s1, a_shu_g1, v_shu_s1, v_shu_g1, store_dadz, store_dvdz, cs=0.2*1e5, time=[1e12, 2e12, 4e12, 8e12]):
        """
        EOM for Terebey, Shu, Cassen (1984)

        Args:
            z_g1, z_s1              : the radius coordinate to calculate upon, in z=1/x
            a_shu_s1                : the Shu (1977) solution for density a0 (array)
            v_shu_s1                : the Shu (1977) solution for velocity v0 (array)
            store_dadz, store_dvdz  : dadz and dvdz from Shu (1977) calculation
            cs                      : sound speed given by S77 to bring dimensionless solution back with dimension      ([cm/s]) 
            time                    : time is provided to bring dimensionless solution back with dimension              ([s])
        """
        self.z_g1       = z_g1
        self.z_s1       = z_s1
        self.a_shu_s1   = a_shu_s1
        self.a_shu_g1   = a_shu_g1
        self.v_shu_s1   = v_shu_s1
        self.v_shu_g1   = v_shu_g1
        self.cs         = cs
        self.time       = time
        
        # quadrupolar EOM for x>1 
        def quadrupolar_greater1_ode(param, z):
            a = param[0]
            v = param[1]
            w = param[2]
            q = param[3]
            p = param[4]

            # equation 69(a)-(e), changing to z=1/x parameter as our boundary condition works for $x\to\infty$
            dadz = 1/(z**2-1) * ((-12*z**2*w) + 2*z**2*(-a/z + 2*v -3*q*z**4 + 2*p/z))
            dvdz = 1/(z**2-1) * (1/z*(-a/z + 2*v - 3*q*z**4 + 2*p/z) - 6*z*w)
            dwdz = -2*w/z -a/2/z**2 + q*z**3 + p/z**2
            dqdz = -5*a/z**6
            dpdz = a/5/z

            # physical constraint, dadz can never be negative, i.e. you can never take away density
            if (dadz < 0):
                dadz = -dadz
            
            return [dadz, dvdz, dwdz, dqdz, dpdz]

        # quadrupolar EOM for x<1
        def quadrupolar_smaller1_ode(z, param, a0, v0, da0dz, dv0dz):
            a = param[0]
            v = param[1]
            w = param[2]
            q = param[3]
            p = param[4]

            m0 = 1/z**2 * a0 * (1/z-v0)
            psi = -1/(3*z**2) - q*z**3 - p/z**2
            dpsidz = 2/(3*z) - 3*q*z**2 + 2*p/z

            A = a*z**2*(2*v0/z - dv0dz) + v*z**2*(2*a0/z - da0dz) - 6*z*a0*w
            B = -a/(a0**2)*da0dz*z**2 + (2-z**2*dv0dz)*v - z**2*dpsidz - 2/3*(m0/2)**4*z**3

            # starting from equation 69, 1. change to z=1/x parameter, 2. adapt boundary conditions from equation 26(b), 67(b), and 68\
            dadz = -1/((1-z*v0)**2 - z**2) * ((1/z-v0)*A + a0*B)
            dvdz = -1/((1-z*v0)**2 - z**2) * ((1/z-v0)*B + A/a0)
            dwdz = 1/(z*v0-1) * ((2/z+v0)*w + a/a0 + psi + (m0/2)**4/3*z**2)
            dqdz = -a/5/z**6
            dpdz = -a/5/z

            # progress recorder
            print("quadruploar x<1, z:", z)

            return [dadz, dvdz, dwdz, dqdz, dpdz]

        # monopolar EOM for x>1
        def monopolar_greater1_ode(param, z): 
            a = param[0]
            v = param[1]
            m = param[2]

            dpsidz = -m - 1/(6*z**3)

            B = a/z + 2*v -z**2*dpsidz - 2/(3*z)

            dadz = 1/(z**2-1) * (a*B)
            dvdz = 1/(z**2-1) * (B/z)
            dmdz = -1/z**4 * (a-1/2)

            return [dadz, dvdz, dmdz]

        # monopolar EOM for x<1
        def monopolar_smaller1_ode(z, param, a0, v0, da0dz, dv0dz):
            a = param[0]
            v = param[1]
            m = param[2]

            m0 = 1/z**2 * a0 * (1/z-v0)
            dpsidz = -m - 1/(6*z**3)

            # 1. change to z=1/x parameter, 2. following equation 73 and 74, following constraints of equation 78, 79
            A = a*z**2*(2*v0/z - dv0dz) + v*z**2*(2*a0/z - da0dz)
            B = a/(a0**2)*z**2*da0dz + (2-z**2*dv0dz)*v - z**2*dpsidz - 2/3*(m0/2)**4*z**3

            dadz = 1/((1-z*v0)**2 - z**2) * ((1/z-v0)*A + a0*B)
            dvdz = -1/((1-z*v0)**2 - z**2) * ((1/z-v0)*B + A/a0)
            dmdz = -1/z**4 * (a-1/2)

            # cheating?    
            #if (dadz < 0):
            #    dvdz = -dadz

            print("monopolar x<1, z:", z)
            
            return [dadz, dvdz, dmdz]

        # boundary condition for quandrupolar, x>1 model
        # adopting equation 65 and 70, K given in III-c of the paper
        K = 1.5*10**(-3)
        # a, v, w, q, p
        param0_g1 = [0, 0, 0, K, 0]

        # solution for quadrupolar, x>1 model
        sol_g1 = odeint(quadrupolar_greater1_ode, param0_g1, z_g1)
        self.a_g = sol_g1[:, 0]
        self.v_g = sol_g1[:, 1]
        self.w_g = sol_g1[:, 2] #continuous through x=1

        # boaundary condition for quadrupolar, x<1 model
        # adopting equation 67 and 68 to go over the x=1 discontinuity. a_shu and v_shu from Shu(1976)
        # here, we manually take limit by sampling points close to 1, such that a(1+e) = a(1) = a(1-e), where e is an arbitrary small number
        # a, v, w, q, p
        param0_s1 = [self.a_g[-1]+2/3*self.v_g[-1], 2/3*self.v_g[-1], self.w_g[-1], K*(1+1/35), -2*K/245]

        # solution for quadrupolar, x<1 model
        sol_s1 = np.zeros((len(z_s1), len(param0_s1)))
        sol_s1[0, :] = param0_s1
        f_s1 = ode(quadrupolar_smaller1_ode)
        f_s1.set_initial_value(param0_s1, z_s1[0])
        f_s1.set_integrator("vode", method = "bdf", atol=1e-20, nsteps=1000)
        f_s1.set_f_params(a_shu_s1[0], v_shu_s1[0], store_dadz[0], store_dvdz[0])

        for i in range(1, len(z_s1)):
            sol_s1[i, :] = f_s1.integrate(z_s1[i])
            f_s1.set_initial_value(sol_s1[i, :], z_s1[i])
            f_s1.set_f_params(a_shu_s1[i], v_shu_s1[i], store_dadz[i], store_dvdz[i])

        self.a_s = sol_s1[:, 0]
        self.v_s = sol_s1[:, 1]
        self.w_s = sol_s1[:, 2]

        param0_mono_g1 = [1/2, 0, 0]
        # solution for monopolar, x>1 model
        sol_mono_g1 = odeint(monopolar_greater1_ode, param0_mono_g1, z_g1)
        self.a_mono_g = sol_mono_g1[:, 0]
        self.v_mono_g = sol_mono_g1[:, 1]
        self.m_mono_g = sol_mono_g1[:, 2] #continuous through x=1

        # initial condition for monopolar x<1 model
        param0_mono_s1 = [self.a_mono_g[-1], self.v_mono_g[-1], self.m_mono_g[-1]]

        # solution for monopolar, x<1 model
        sol_mono_s1 = np.zeros((len(z_s1), len(param0_mono_s1)))
        sol_mono_s1[0, :] = param0_mono_s1
        f_mono_s1 = ode(monopolar_smaller1_ode)
        f_mono_s1.set_integrator("vode", method = "bdf", atol=1e-20, nsteps=1000)
        f_mono_s1.set_initial_value(param0_mono_s1, z_s1[0])
        f_mono_s1.set_f_params(a_shu_s1[0], v_shu_s1[0], store_dadz[0], store_dvdz[0])

        for i in range(1, len(z_s1)):
            sol_mono_s1[i, :] = f_mono_s1.integrate(z_s1[i])
            f_mono_s1.set_initial_value(sol_mono_s1[i, :], z_s1[i])
            f_mono_s1.set_f_params(a_shu_s1[i], v_shu_s1[i], store_dadz[i], store_dvdz[i])

        self.a_mono_s1 = [i for i in sol_mono_s1[:, 0]]
        self.v_mono_s1 = [i for i in sol_mono_s1[:, 1]]

    def get_r_g1(self):
        r = []
        for t in self.time:
            r.append([1/i*self.cs*t for i in self.z_g1])
        
        return r

    def get_r_s1(self):
        r = []
        for t in self.time:
            r.append([1/i*self.cs*t for i in self.z_s1])
        
        return r

    ##### Density #####
    def get_quad_a_g1(self):
        return self.a_g
    
    def get_quad_a_s1(self):
        return self.a_s
    
    def get_mono_a_g1(self):
        return self.a_mono_g
    
    def get_mono_a_s1(self):
        return self.a_mono_s1
    
    ##### Velocity #####
    def get_quad_v_g1(self):
        return self.v_g
    
    def get_quad_v_s1(self):
        return self.v_s
    
    def get_quad_w_g1(self):
        return self.w_g
    
    def get_quad_w_s1(self):
        return self.w_s
    
    def get_mono_v_g1(self):
        return self.v_mono_g
    
    def get_mono_v_s1(self):
        return self.v_mono_s1
    
    def get_cgs_v_r_g1(self, tau=1e-3):
        leg = legendre(2)
        
        v_2 = []
        theta_range = np.linspace(0, 2*np.pi, param.theta_pts)
        for theta in range(len(theta_range)):
            v_2_theta = []
            for i in range(len(self.v_g)):
                v_2_theta_i = self.v_mono_g[i] + self.v_g[i]*leg(np.cos(theta))
                v_2_theta.append(v_2_theta_i)
            v_2.append(v_2_theta)
            
        v_r = []
        for i in range(len(v_2)):
            v_j = []
            for j in range(len(v_2[i])):
                v_i = self.v_shu_g1[j] + tau*v_2[i][j]
                v_j.append(v_i)
            v_r.append(v_j)

        vel_rad = []
        for t in self.time:
            vel_rad.append([])
            for i in range(len(v_r)):
                vel_rad[-1].append([self.cs*j for j in v_r[i]])
        
        return vel_rad
    
    def get_cgs_v_r_s1(self, tau=1e-3):
        leg = legendre(2)
        
        v_2 = []
        theta_range = np.linspace(0, 2*np.pi, param.theta_pts)
        for theta in range(len(theta_range)):
            v_2_theta = []
            for i in range(len(self.v_g)):
                v_2_theta_i = self.v_mono_s1[i] + self.v_s[i]*leg(np.cos(theta))
                v_2_theta.append(v_2_theta_i)
            v_2.append(v_2_theta)
            
        v_r = []
        for i in range(len(v_2)):
            v_j = []
            for j in range(len(v_2[i])):
                v_i = self.v_shu_s1[j] + tau*v_2[i][j]
                v_j.append(v_i)
            v_r.append(v_j)

        vel_rad = []
        for t in self.time:
            vel_rad.append([])
            for i in range(len(v_r)):
                vel_rad[-1].append([self.cs*j for j in v_r[i]])
        
        return vel_rad

    def get_cgs_v_theta_g1(self, tau=1e-3):
        leg = legendre(2)
        leg_poly = Polynomial(leg)
        dev_leg_poly = leg_poly.deriv()

        w_2 = []
        theta_range = np.linspace(0, 2*np.pi, param.theta_pts)
        for theta in range(len(theta_range)):
            w_2_theta = []
            for i in range(len(self.w_g)):
                w_2_theta_i = self.w_g[i] * (-dev_leg_poly(np.cos(theta))*np.sin(theta))
                w_2_theta.append(w_2_theta_i)
            w_2.append(w_2_theta)
        
        v_theta = []
        for i in range(len(w_2)):
            w_j = []
            for j in range(len(w_2[i])):
                w_i = tau*w_2[i][j]
                w_j.append(w_i)
            v_theta.append(w_j)

        vel_theta = []
        for t in self.time:
            vel_theta.append([])
            for i in range(len(v_theta)):
                vel_theta[-1].append([self.cs*j for j in v_theta[i]])
        
        return vel_theta
    
    def get_cgs_v_theta_s1(self, tau=1e-3):
        leg = legendre(2)
        leg_poly = Polynomial(leg)
        dev_leg_poly = leg_poly.deriv()

        w_2 = []
        theta_range = np.linspace(0, 2*np.pi, param.theta_pts)
        for theta in range(len(theta_range)):
            w_2_theta = []
            for i in range(len(self.w_s)):
                w_2_theta_i = self.w_s[i] * (-dev_leg_poly(np.cos(theta))*np.sin(theta))
                w_2_theta.append(w_2_theta_i)
            w_2.append(w_2_theta)
        
        v_theta = []
        for i in range(len(w_2)):
            w_j = []
            for j in range(len(w_2[i])):
                w_i = tau*w_2[i][j]
                w_j.append(w_i)
            v_theta.append(w_j)

        vel_theta = []
        for t in self.time:
            vel_theta.append([])
            for i in range(len(v_theta)):
                vel_theta[-1].append([self.cs*j for j in v_theta[i]])
        
        return vel_theta
    
    def get_cgs_v_phi(self, tau=1e-3):
        x   = [1/i for i in self.z_s1][::-1]+[1/i for i in self.z_g1]
        a0  = [i for i in self.a_shu_s1][::-1]+[i for i in self.a_shu_g1][::-1]
        v0  = [i for i in self.v_shu_s1][::-1]+[i for i in self.v_shu_g1][::-1]

        m0 = []
        for i in range(len(x)):
            m0.append(x[i]**2 * a0[i] * (x[i]-v0[i]))

        sigma_2 = []
        theta_range = np.linspace(0, 2*np.pi, param.theta_pts)
        for theta in range(len(theta_range)):
            sigma_2_theta = []
            for i in range(len(m0)):
                sigma_2_theta_i = (1/2*m0[i] * np.sin(theta))**2
                sigma_2_theta.append(sigma_2_theta_i)
            sigma_2.append(sigma_2_theta)

        sigma = []
        for i in range(len(sigma_2)):
            sigma_j = []
            for j in range(len(sigma_2[i])):
                sigma_j_i = tau * sigma_2[i][j]
                sigma_j.append(sigma_j_i)
            sigma.append(sigma_j)
        
        radius = []
        for i in range(len(self.time)):
            radius.append(self.get_r_s1()[i][::-1]+self.get_r_g1()[i][::-1])

        vel_phi = []
        for t in range(len(self.time)):
            vel_phi.append([])
            for theta in range(len(theta_range)):
                vel_phi_theta = []
                for i in range(len(radius)):
                    vel_phi_i = self.cs**2/param.omega/(radius[t][i]*np.sin(theta)) * sigma[theta][i]
                    vel_phi_theta.append(vel_phi_i)
                vel_phi[-1].append(vel_phi_theta)
        
        return vel_phi

    def get_cgs_density_g1(self, tau=1e-3):
        # [time][theta]
        # for omega=1e-14 -> tau^2=1e-3
        a_2 = []
        leg = legendre(2)

        theta_range = np.linspace(0, 2*np.pi, param.theta_pts)
        for theta in range(len(theta_range)):
            a_2_theta = []
            for i in range(len(self.a_g)):
                a_2_theta_i = self.a_mono_g[i] + self.a_g[i]*leg(np.cos(theta))
                a_2_theta.append(a_2_theta_i)
            a_2.append(a_2_theta)
        
        rho = []
        for i in range(len(a_2)):
            rho_j = []
            for j in range(len(a_2[i])):
                rho_i = self.a_shu_g1[j] + tau*a_2[i][j]
                rho_j.append(rho_i)
            rho.append(rho_j)

        density = []
        for t in self.time:
            density.append([])
            for i in range(len(rho)):
                density[-1].append([(j/(4*np.pi*param.G*t**2)) for j in rho[i]])

        return density
    
    def get_cgs_density_s1(self, tau=1e-3):
        # [time][theta]
        # for omega=1e-14 -> tau^2=1e-3
        a_2 = []
        leg = legendre(2)

        theta_range = np.linspace(0, 2*np.pi, param.theta_pts)
        for theta in range(len(theta_range)):
            a_2_theta = []
            for i in range(len(self.a_s)):
                a_2_theta_i = self.a_mono_s1[i] + self.a_s[i]*leg(np.cos(theta))
                a_2_theta.append(a_2_theta_i)
            a_2.append(a_2_theta)
        
        rho = []
        for i in range(len(a_2)):
            rho_j = []
            for j in range(len(a_2[i])):
                rho_i = self.a_shu_s1[j] + tau*a_2[i][j]
                rho_j.append(rho_i)
            rho.append(rho_j)

        density = []
        for t in self.time:
            density.append([])
            for i in range(len(rho)):
                density[-1].append([(j/(4*np.pi*param.G*t**2)) for j in rho[i]])

        return density
    
    def get_nh_density_g1(self, tau=1e-3):
        # [time][theta]
        # for omega=1e-14 -> tau^2=1e-3
        a_2 = []
        leg = legendre(2)

        theta_range = np.linspace(0, 2*np.pi, param.theta_pts)
        for theta in range(len(theta_range)):
            a_2_theta = []
            for i in range(len(self.a_g)):
                a_2_theta_i = self.a_mono_g[i] + self.a_g[i]*leg(np.cos(theta))
                a_2_theta.append(a_2_theta_i)
            a_2.append(a_2_theta)
        
        rho = []
        for i in range(len(a_2)):
            rho_j = []
            for j in range(len(a_2[i])):
                rho_i = self.a_shu_g1[j] + tau*a_2[i][j]
                rho_j.append(rho_i)
            rho.append(rho_j)

        density = []
        for t in self.time:
            density.append([])
            for i in range(len(rho)):
                density[-1].append([(j/(4*np.pi*param.G*t**2))/param.M_H for j in rho[i]])

        return density
    
    def get_nh_density_s1(self, tau=1e-3):
        # [time][theta]
        # for omega=1e-14 -> tau^2=1e-3
        a_2 = []
        leg = legendre(2)

        theta_range = np.linspace(0, 2*np.pi, param.theta_pts)
        for theta in range(len(theta_range)):
            a_2_theta = []
            for i in range(len(self.a_s)):
                a_2_theta_i = self.a_mono_s1[i] + self.a_s[i]*leg(np.cos(theta))
                a_2_theta.append(a_2_theta_i)
            a_2.append(a_2_theta)
        
        rho = []
        for i in range(len(a_2)):
            rho_j = []
            for j in range(len(a_2[i])):
                rho_i = self.a_shu_s1[j] + tau*a_2[i][j]
                rho_j.append(rho_i)
            rho.append(rho_j)

        density = []
        for t in self.time:
            density.append([])
            for i in range(len(rho)):
                density[-1].append([(j/(4*np.pi*param.G*t**2))/param.M_H for j in rho[i]])

        return density