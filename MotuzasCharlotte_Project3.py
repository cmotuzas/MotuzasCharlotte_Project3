# Charlotte Motuzas 
# Project 3 

import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as scint  
import csv
import astropy.units as astro 

# Question 1 

# Use the Python ODE integrator solve_ivp to solve Equations (8) and (9) subject to the boundary conditions 
# m(r=0) = 0 and ρ(r=0) = ρc. Your goal is to determine the radius and mass of model white dwarfs by starting 
# the integration at r = 0 and integrating outward to a distance rs where the density drops to zero, 
# i.e. , ρ(rs)= 0. Calculate solutions for 10 values of ρc in the range 10−1 to 2.5 · 106 and take μe=2. 
# [Hint: if starting your integration at exactly r=0 raises a ZeroDivisionError, can you start it at some tiny 
# value of r instead? How do you know what 'tiny' is in this context?

# use events function to find when density is zero
# 1e-5 or something 

def dstatedt(r,statevec):
    rho = statevec[0]
    m = statevec[1]
    x = np.cbrt(rho)
    gamma = (x**2)/(3*np.sqrt((1+x**2)))
    return [-m*rho/(gamma*(r**2)),rho*(r**2)]

def stop_event(r,statevec):
    if statevec[0] - 1e-5 < 0:
        statevec[0] = 0
    return statevec[0]

stop_event.terminal = True


rho_c = np.logspace(-1,6.3,10)
mu_e = 2
M0 = (5.67e33/(mu_e**2))/1000
R0 = 7.72e8/mu_e
r_i = 1e-10 # told to start at a tiny version of r
rho_0 = (9.74*1e5)*mu_e # g/cm^3
rs = 10
r_soln = np.linspace(r_i,rs,10000)
M = np.empty(len(rho_c))
R = np.empty(len(rho_c))
M_msun = np.empty(len(rho_c))
R_rsun = np.empty(len(rho_c))
for i in range(len(rho_c)):
    v0 = np.array([rho_c[i],0])
    result = scint.solve_ivp(dstatedt,(r_i,rs),v0,t_eval=r_soln,events=stop_event)
    
    # Question 2 - Transforming into physical units
    M[i] = result.y_events[0][0][1]*M0 # kg 
    R[i] = result.t_events[0][0]*R0*1e-5
    M_msun[i] = result.y_events[0][0][1]*M0*4/1.989e30 # kg 
    R_rsun[i] = result.t_events[0][0]*R0*1e-5/696340

    plt.scatter(M[i],R[i],label="$\\rho_c$ = {}".format(round(rho_c[i],1)))
plt.title('White Dwarf Mass vs Radius')
plt.xlabel('Mass (kg)')
plt.ylabel('Radius (km)')
plt.legend()
plt.show()

for i in range(len(rho_c)): # done as a loop for labelling purposes
    plt.scatter(M_msun[i],R_rsun[i],label="$\\rho_c$ = {}".format(round(rho_c[i],1)))
plt.title('White Dwarf Mass vs Radius')
plt.xlabel('Mass ($M_{sun}$)')
plt.ylabel('Radius (km)')
#plt.legend()

# Question 3 

# Pick 3 values of ρc and run solve_ivpagain, choosing a different integration method (of the ones available in solve_ivp) 
# from the one you used in part (1). How different are your results? 

# RK23

M_23 = np.empty(3)
R_23 = np.empty(3)
diff = np.empty(3)
for i in range(3):
    v0 = np.array([rho_c[3*i],0])
    result = scint.solve_ivp(dstatedt,(r_i,rs),v0,method='RK23',t_eval=r_soln,events=stop_event)
    
    # Question 2 - Transforming into physical units
    M_23[i] = result.y_events[0][0][1]*M0 # kg 
    R_23[i] = result.t_events[0][0]*R0*1e-5

    diff[i] = np.abs((M[3*i]-M_23[i])*100/M[3*i])

print(diff)

# Question 4 

#Tremblay et al. (2017) give white dwarf masses and radii obtained from binary star systems with distances measured by the 
# Gaia satellite. These data are in the file wd_mass_radius.csv, available on here (measurements are in units of the Sun's 
# mass and radius). Plot these observed data with their error bars on your computed mass-radius relation, paying attention 
# to the units. How well do the observations agree with your calculations?

M_Msun = []
M_unc = []
R_Rsun = []
R_unc = []

with open('wd_mass_radius.csv') as csvfile: 
    csvReader = csv.reader(csvfile)
    next(csvReader, None)  # skip the headers
    for row in csvReader: 
        M_Msun.append(float(row[0]))
        M_unc.append(float(row[1]))
        R_Rsun.append(float(row[2]))
        R_unc.append(float(row[3]))

print(M_Msun)
print(R_Rsun)


plt.errorbar(M_Msun,R_Rsun,yerr=R_unc,xerr=M_unc,fmt='o',linewidth=0.5)
plt.show()

#plt.errorbar(np.array([10,40,70]),np.array([np.average(T_exp_10),np.average(T_exp_40),np.average(T_exp_70)]),yerr=errorvec,fmt='-',linewidth=0.5)

