# Charlotte Motuzas 
# Project 3 

import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as scint  
import csv
import astropy.units as astro 
from decimal import Decimal


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
    '''This function returns the derivaives of the state variables rho (density, dimensionless equivalent) and m (mass, dimensionless equivalent). 
    This function takes inputs r and statevec, where statevec is a two dimensional vector referring to rho and mass respectively. This function is 
    to be utilized within the solve_ivp() function.'''
    rho = statevec[0]
    m = statevec[1]
    x = np.cbrt(rho)
    gamma = (x**2)/(3*np.sqrt((1+x**2)))
    return [-m*rho/(gamma*(r**2)),rho*(r**2)]

def stop_event(r,statevec):
    '''This function is used within the solve_ivp() function to identify and terminate when the solution returns a density less than 1e-5.
    This takes input variables r and statevec, where r is the independant variable and statevec is a two dimensional vector referring to rho and mass respectively.
    This function makes any value of density equal to zero if less than 1e-5, and then returns the density value to be identified within the 'events' condition in solve_ivp()'''
    if statevec[0] - 1e-5 < 0:
        statevec[0] = 0
    return statevec[0]

stop_event.terminal = True # event stopping function terminates the solution 


rho_c = np.logspace(-1,6.3,10) # range of central densities in logspace 
mu_e = 2
M0 = (5.67e33/(mu_e**2))/1000 # conversion constant for mass
R0 = 7.72e8/mu_e # conversion constant for radius 
r_i = 1e-10 # told to start at a tiny version of r
rho_0 = (9.74*1e5)*mu_e # g/cm^3, conversion constant for density 
rs = 10 
r_soln = np.linspace(r_i,rs,10000) # values for which the solution should be computed
Mch = (5.836/(mu_e**2))*1.989e30
# initializing mass and radius arrays 
M = np.empty(len(rho_c))
R = np.empty(len(rho_c))
for i in range(len(rho_c)):
    # looping through each of the center density values to obtain ten solutions for internal density structure 
    v0 = np.array([rho_c[i],0])
    result = scint.solve_ivp(dstatedt,(r_i,rs),v0,t_eval=r_soln,events=stop_event)
    
    # Question 2 - Transforming into physical units, plotting, comparing to Mch
    M[i] = result.y_events[0][0][1]*M0 # kg 
    R[i] = result.t_events[0][0]*R0*1e-5
                                            
    plt.scatter(M[i],R[i],label=f'$\\rho_c$ = {Decimal(rho_c[i]*rho_0):.2e} $g/cm^3$')
plt.plot(np.array([Mch,Mch]),np.array([-200,15000]),'k--',label="Chandrasekhar Limit",linewidth=0.8)
plt.title('White Dwarf Mass vs Radius')
plt.xlabel('Mass (kg)')
plt.ylabel('Radius (km)')
plt.ylim([-200,15000])
plt.legend(fontsize=9)
plt.show()

# Question 3 

# Pick 3 values of ρc and run solve_ivpagain, choosing a different integration method (of the ones available in solve_ivp) 
# from the one you used in part (1). How different are your results? 

# RK23

# initializing arrays for mass and radius for the second method 
M_23 = np.empty(3)
R_23 = np.empty(3)
diff = np.empty(3)
diffR = np.empty(3)
for i in range(3):
    v0 = np.array([rho_c[3*i],0])
    result = scint.solve_ivp(dstatedt,(r_i,rs),v0,method='RK23',t_eval=r_soln,events=stop_event) # third-order runge-kutta method 
    
    # Question 2 - Transforming into physical units
    M_23[i] = result.y_events[0][0][1]*M0 # kg 
    R_23[i] = result.t_events[0][0]*R0*1e-5

    diff[i] = np.abs((M[3*i]-M_23[i])*100/M[3*i])
    diffR[i] = np.abs((R[3*i]-R_23[i])*100/R[3*i])


print('For rho_c {}, the second method produces a mass that varies by {}, and a radius that differs by {}'.format(round(rho_c[0],2),round(diff[0],2),round(diffR[0],2)))
print('For rho_c {}, the second method produces a mass that varies by {}, and a radius that differs by {}'.format(round(rho_c[3],2),round(diff[1],2),round(diffR[1],2)))
print('For rho_c {}, the second method produces a mass that varies by {}, and a radius that differs by {}'.format(round(rho_c[6],2),round(diff[2],2),round(diffR[2],2)))

# Question 4 

#Tremblay et al. (2017) give white dwarf masses and radii obtained from binary star systems with distances measured by the 
# Gaia satellite. These data are in the file wd_mass_radius.csv, available on here (measurements are in units of the Sun's 
# mass and radius). Plot these observed data with their error bars on your computed mass-radius relation, paying attention 
# to the units. How well do the observations agree with your calculations?

rho_c = np.logspace(-1,6.3,50) # larger central density vector for initial conditions

# initializing arrays for mass and radius 
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
    M_msun[i] = result.y_events[0][0][1]*M0/1.989e30 # kg 
    R_rsun[i] = result.t_events[0][0]*R0*1e-5/696340


# initializing the arrays prior to importing the CSV file of Gaia data 
M_Msun_data = []
M_unc = []
R_Rsun_data = []
R_unc = []

with open('wd_mass_radius.csv') as csvfile: 
    csvReader = csv.reader(csvfile) # reading the file 
    next(csvReader, None)  # skip the headers
    for row in csvReader: 
        # assigning csv data to vectors 
        M_Msun_data.append(float(row[0]))
        M_unc.append(float(row[1]))
        R_Rsun_data.append(float(row[2]))
        R_unc.append(float(row[3]))

plt.plot(M_msun,R_rsun,'r')
plt.errorbar(M_Msun_data,R_Rsun_data,yerr=R_unc,xerr=M_unc,fmt='o',linewidth=0.5,markersize=4)
plt.title('White Dwarf Mass vs Radius')
plt.xlabel('Mass ($M_{sun}$)')
plt.ylabel('Radius ($R_{sun}$)')
plt.legend(['Computed Mass-Radius Relation','Gaia Measurements'])
plt.show()

