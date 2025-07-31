# I have no idea what is happening. Probably useless now


import numpy as np
import matplotlib.pyplot as plt
import math
from Functions import dr_dlambda, dphi_dlambda, L_init,pt_init
from Functions import calculate_r_1,calculate_r_2,calculate_r_3
from scipy.integrate import solve_ivp

# Parameters
M, L= 1.0, 2.0
theta=math.pi/2
B=0.10

# Initial values of r and L (calculated from null condition)
r_1=calculate_r_1(B,M,theta)
r_3=calculate_r_3(B,M,theta)
# L=L=L_init(r_1,theta,B,M,E)
E=pt_init(L,theta,M,B,r_3)
# Range of affine parameter
lambda_span = (0, 20)  

# Solve the ODE
sol = solve_ivp(dr_dlambda, lambda_span, [r_3], args=(E,B,L,M,1), method='RK45', t_eval=np.linspace(0, 20, 100))
sol2 = solve_ivp(dphi_dlambda, lambda_span, [r_1], args=(B,L), method='RK45', t_eval=np.linspace(0, 20, 100))

# Plots r vs lambda
plt.plot(sol.t, sol.y[0], label="Radial evolution")
plt.xlabel(r"Affine parameter $\lambda$")
plt.ylabel("$r$")
plt.title("Evolution of $r$ from stable turning point $r_1$")
plt.legend()
plt.show()

# Plot phi vs r
# plt.plot(sol.y[0], sol2.y[0], label="Radial evolution")
# plt.ylabel("$\phi$")
# plt.xlabel("$r$")
# plt.title("Evolution of $r$ and $\phi$ from stable turning point $r_1$")
# plt.legend()
# plt.show()