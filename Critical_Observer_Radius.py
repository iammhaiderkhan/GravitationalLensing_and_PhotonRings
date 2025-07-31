# Calculates and plots the critical observer radius
# beyond which the horizontal shadow component formulas
# is invalid for an equaotorial observer


import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import brentq
from Functions import impact_parameter, calculate_r_3, calculate_r_1, Lambda

theta=math.pi/2
M = 1

def RHS(r,B):
    term1 = Lambda(r,theta,B)**4*(r-2*M)/r**3
    r2=calculate_r_3(B,M,theta)
    term2 = r2**3/((Lambda(r2,theta,B)**4*(r2-2*M)))
    return term1*term2

def sin2beta(r_O, r_2, M, B):
    """Compute sin^2(beta) for horizontal shadow radius."""
    lam_O = Lambda(r_O,theta, B)
    lam_2 = Lambda(r_2,theta, B)
    numerator = lam_O**4 * (r_O - 2*M) * r_2**3
    denominator = r_O**3 * lam_2**4 * (r_2 - 2*M)
    return numerator / denominator

def find_r_O_prime(B, M):
    """
    For a given B, solve sin^2(beta) = 1 to find r_O'.
    - B: magnetic field parameter
    - M: black hole mass
    - calculate_r_2: function that returns the unstable photon radius r_2 for given B, M
    """
    r_2 = calculate_r_3(B,M,theta)
    # Define the equation sin^2(beta) - 1 = 0
    def f(r_O):
        return sin2beta(r_O, r_2, M, B) - 1
    
    # Bracket for root finding: start just above r_2 and extend to large radius
    r_min = r_2 * 1.01
    r_max = 1000 * M
    
    # Use Brent's method to find r_O'
    r_O_prime = brentq(f, r_min, r_max)
    print(f"r_O' for B = {B:.2f}M is {r_O_prime:.3f}M")
    print(f"sin^2 beta = {RHS(r_O_prime,B):.2f}")
    
    return r_O_prime

def plot_r_O_vs_B(B_vals, M):
    """
    Plot r_O' as a function of B.
    - B_vals: array of B values
    - M: black hole mass
    - calculate_r_2: function to compute unstable photon radius r_2
    """
    r_O_primes = []
    for B in B_vals:
        r_O_primes.append(find_r_O_prime(B, M))
    
    plt.figure(figsize=(8, 6))
    plt.plot(B_vals, r_O_primes)
    plt.axhline(calculate_r_1(0.189,M,theta), color="gray", linestyle="--", linewidth=1, label=fr"$r_1=r_2=$ {calculate_r_3(0.189,M,theta):.3f}M for $B_c=0.189/M$")
    plt.xlabel("Magnetic Field Strength $B$",fontsize=12)
    plt.ylabel(r"Critical Observer Radius $r_O'$",fontsize=12)
    # plt.title(r"Maximum Observer Radius $r_O'$ vs Magnetic Field $B$",fontsize=14)
    plt.legend()
    plt.grid(True)
    plt.show()

# Example usage:
# from your_module import calculate_r_2
B_values = np.linspace(0.01, 0.189, 40)
plot_r_O_vs_B(B_values, M=1.0)
# B_c=0.001
# print(f"r_1 for B = {B_c}M = {calculate_r_1(B_c,M,theta):.3f}")
# find_r_O_prime(B_c,M)

# Uncomment and adjust the above lines according to your calculate_r_2 implementation.

B=[0.05,0.1,0.15]
r=[calculate_r_3(i,M,theta) for i in B]
print(r)