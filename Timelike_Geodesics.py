import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math
from Functions import Lambda
from Functions import calculate_r_1, calculate_r_3, theta_dot_squared2, cubic_equation


def expression(r,B,M):
    return 1/4*B**4*r**5 + 4*M*r**2*B**2 - 4*r + 8*M + M*B**4*r**4

def E_init(r,theta,B,M):
    Lmbd=Lambda(r,theta,B)
    num = 2 * r - (r**3 * B**2 / Lmbd) - (2 * M * r**2 * B**2 / Lmbd) - (4 * M / Lmbd**2)
    den = (2 * r / Lmbd**2) - (2 * r**3 * B**2 / Lmbd**3) - (2 * M / Lmbd**2)/(1-2*M/r)
    return (num/den)
    # return math.sqrt(num/den)

def L_init(r,B,M):
    term1= r**2/Lambda(r,math.pi/2,B)**3
    term3= expression(r,B,M)
    term4=cubic_equation(r,B,M,math.pi/2)*Lambda(r,math.pi/2,B)
    term5=r**2/Lambda(r,math.pi/2,B)**2
    return term1*term3/term4-term5

def Double_derivative(r,M,B,L,E):
    numerator = (
        -256 * (
        L**2 * (B**2 * r**2 + 4)**4 * (
            7 * B**4 * r**3 * (8 * M - 3 * r)
            + 24 * B**2 * r * (3 * r - 4 * M)
            - 16
        )
        + 96 * r * (
            M * (15 * B**4 * r**4 - 48 * B**2 * r**2 + 16) * (B**2 * r**2 + 4)**2
            - 2 * r * (3 * B**4 * r**4 - 14 * B**2 * r**2 + 8) * (B**2 * r**2 + 4)**2
            + 16 * E**2 * r * (13 * B**4 * r**4 - 40 * B**2 * r**2 + 16)
        )
    ))

    denominator = L**2 * (B**2 * r**2 + 4)**10

    return numerator / denominator

M = 1.0
B = 0.15
theta = math.pi/2

r1=calculate_r_1(B,M,theta)
r2=calculate_r_3(B,M,theta)

# print(L_init(r1+0.01,B,M))
# print(L_init(r1-0.01,B,M))


# Define the range and resolution
r_values = np.linspace(3, 30, 1000)
eps = 1e-2     # Small number to exclude values too close to r1

# Filter out values too close to r1
# r_filtered = r_values[
#     (np.abs(r_values - r1) > eps) 
#   & (np.abs(r_values - r2) > eps)
# ]


# Evaluate L_init at safe r-values
L_values = []
for r in r_values:
    try:
        L = L_init(r, B, M)
        if np.isfinite(L):  # Avoid inf/nan due to numerical issues
            L_values.append(L)
        else:
            L_values.append(np.nan)
    except:
        L_values.append(np.nan)

# Plotting
plt.figure(figsize=(8, 6))
plt.plot(r_values, L_values, color="blue")

plt.axvline(r2, color="#d62728", linestyle="--", linewidth=2.5, label=rf"$r_2 = {r2:.3f}M$")
plt.axvline(r1, color="#9467bd", linestyle="--", linewidth=2.5, label=rf"$r_1 = {r1:.3f}M$")
plt.axhline(0, color="gray", linestyle="--", linewidth=1)

plt.xlabel("r",fontsize=12)
plt.ylabel("$L^2$",fontsize=12)
plt.ylim(-500,500)
# plt.title("Plot of $L^2$ excluding the singularity at $r_1$")
plt.legend()
plt.grid(True)
plt.show()
