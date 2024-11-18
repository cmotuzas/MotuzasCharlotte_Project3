# Charlotte Motuzas 
# Project 3 

import numpy as np 
import matplotlib.pyplot as plt 
import scipy.integrate as scint  
import scipy.interpolate as interp
import astropy.units as astro 

# Question 1 

# Use the Python ODE integrator solve_ivp to solve Equations (8) and (9) subject to the boundary conditions 
# m(r=0) = 0 and ρ(r=0) = ρc. Your goal is to determine the radius and mass of model white dwarfs by starting 
# the integration at r = 0 and integrating outward to a distance rs where the density drops to zero, 
# i.e. , ρ(rs)= 0. Calculate solutions for 10 values of ρc in the range 10−1 to 2.5 · 106 and take μe=2. 
# [Hint: if starting your integration at exactly r=0 raises a ZeroDivisionError, can you start it at some tiny 
# value of r instead? How do you know what 'tiny' is in this context?

# Question 2

# Transform your results from the ODE solution into physical (not dimensionless) units and plot R as a function of 
# M(i.e. R on the y-axis and M on the x-axis). Can you estimate the Chandrasekhar limit (the upper limit to white 
# dwarf masses) from your plot? How does your result compare with Kippenhahn & Weigert (1990) who cite MCh = 5.836/(μe)2 ? 
# [Hint: you may find the astropy.units package useful for doing the unit conversions.]

# Question 3 

# Pick 3 values of ρc and run solve_ivpagain, choosing a different integration method (of the ones available in solve_ivp) 
# from the one you used in part (1). How different are your results? 

# Question 4 

#Tremblay et al. (2017) give white dwarf masses and radii obtained from binary star systems with distances measured by the 
# Gaia satellite. These data are in the file wd_mass_radius.csv, available on here (measurements are in units of the Sun's 
# mass and radius). Plot these observed data with their error bars on your computed mass-radius relation, paying attention 
# to the units. How well do the observations agree with your calculations?