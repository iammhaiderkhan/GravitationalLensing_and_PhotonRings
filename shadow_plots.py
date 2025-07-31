# PLots the horizontal and vertical components of the
# shadow as seen by an observer at radius r_O

import numpy as np
import math
import matplotlib.pyplot as plt
from Functions import calculate_r_1, calculate_r_3, impact_parameter

# === 1. Parameters ===
M = 1.0
B = 0.05

# === 3. Helper: Lambda(r) ===
def Lambda(r):
    return 1.0 + 0.25 * B**2 * r**2

r2=calculate_r_3(B,M,math.pi/2)
r1=calculate_r_1(B,M,math.pi/2)

Lam2 = Lambda(r2)

# === 4. Observer radii and angle arrays ===
rO = np.linspace(3.01*M, 50*M, 200)  # from just above horizon out to 200M

# vertical: sin^2 alpha = 27 M^2 (rO-2M)/rO^3
sin2_alpha = 27*M**2 * (rO - 2*M) / rO**3
alpha = np.arcsin(np.sqrt(sin2_alpha))   # radians

# horizontal: sin^2 beta = [Lambda_O^4 (rO-2M)/rO^3] * [r2^3/(Lam2^4 (r2-2M))]
sin2_beta = (Lambda(rO)**4 * (rO - 2*M) / rO**3) * impact_parameter(M,r2,B)**2
beta = np.arcsin(np.sqrt(sin2_beta))

# If i want to see the degree conversions
alpha_deg = np.degrees(alpha)
beta_deg  = np.degrees(beta)


# === 5. Plot ===
plt.figure(figsize=(8,5))
plt.plot(rO/M, alpha, label=r'$\alpha$ (vertical)')
plt.plot(rO/M, beta,  label=r'$\beta$ (horizontal)')
plt.axvline(r2/M, color='gray', linestyle='--', label=r'$r_2$')
plt.axvline(r1/M, color='lightgray', linestyle='--', label=r'$r_1$')

plt.xlabel(r'Observer radius $r_O/M$',fontsize=12)
plt.ylabel('Angle (rad)',fontsize=12)
plt.title(fr'$B = ${B}$/M$',fontsize=14)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
