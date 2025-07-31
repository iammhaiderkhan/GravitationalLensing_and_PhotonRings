# Plots the deflection angle for the meridional plane in the Ernst spacetime
# The formula obtained is the same as in Schwarzschild

from Functions import integrand, deflection_angle
import matplotlib.pyplot as plt
import numpy as np

M = 1.0         
B = 0.1         
E = 1.0     
r_m=np.linspace(2.1,100,500)
delta=[]
for r in r_m:
    delta.append(deflection_angle(r,M,B))

plt.plot(r_m,delta)
plt.xlabel(r"$r_m$")
plt.ylabel(r"$\delta$")
# plt.title(r"")
plt.show()