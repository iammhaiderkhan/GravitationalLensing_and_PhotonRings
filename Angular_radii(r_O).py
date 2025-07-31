# Use case of the impact parameter determination
# and calculation of the shadow angular radii for various observer radii

from Functions import impact_parameter, calculate_r_1, calculate_r_3, deflection_angle, integrand, Lambda
import math
import numpy as np

M = 1.0  # Black hole mass
B = 0.1 # Magnetic Field Strenght
theta = math.pi/2 # Equatorial plane

r_1=calculate_r_1(B,M,theta)
r_3=calculate_r_3(B,M,theta)

b_1=impact_parameter(M,r_1,B)
b_3=impact_parameter(M,r_3,B)

# print(f"Impact parameter for r_1 = {b_1:.4f}")
# print(f"Impact parameter for r_3 = {b_3:.4f}")

# for r_m in np.linspace(2.1, 15, 20):
#     # delta = deflection_angle(r_m, M, B)
#     integrand(10,r_m, M, B)

#     # print(f"r_m = {r_m:.2f}, δ = {delta:.6f} radians")


def expression(r, r_m, M, B):
    Lambda_m = Lambda(r_m,math.pi/2, B)
    Lambda_r = Lambda(r, math.pi/2, B)
    
    # Compute the expression under the square root
    numerator = (Lambda_m**4) / (Lambda_r**8) * (1 - 2 * M / r_m) * r**4
    term2 = r**2 * r_m**2 / Lambda_r**4
    term3 = 2 * M * r * r_m**2 / Lambda_r**4
    
    return numerator - term2 + term3

# for r_m in np.linspace(2.1, 1602.1, 20):
#     print(f"For r= 200 and r_m= {r_m:.4f} = {expression(2000, r_m,M,B)}")

# print(rf"Impact parameter b_3^2 for B = {B}M is {b_3**2:.4f}")

def alpha(r_O,M):
    return 27 * M**2 * (r_O - 2 * M)/r_O**3
def beta(r_O,M,B,b):
    return Lambda(r_O,math.pi/2,B)**4 * b**2 * (r_O - 2 * M)/r_O**3

r_O=[r_3,4,5,6,7,8]
for k in r_O:
    print(f"For r_O = {k}M, sin^2(alpha) = {alpha(k,1)} and sin^2(beta) = {beta(k,1,B,b_3)}")