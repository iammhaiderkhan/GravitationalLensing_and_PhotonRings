import math
import numpy as np
import matplotlib.pyplot as plt
from Functions import Lambda, calculate_r_1, calculate_r_3

M = 1
B_list = [0.01, 0.05, 0.10, 0.18 * M]
theta = math.pi / 2
pert = 0.00001

# Create r_range once for B=0.01
B_base = 0.01
r3_base = calculate_r_3(B_base, M, theta)
r1_base = calculate_r_1(B_base, M, theta)
r_range = np.linspace(r3_base + pert, r1_base - pert, 200)

import numpy as np

def potential(r, theta, M, B, E, L):
    """
    Evaluate -2 * V_{E,L} in the Ernst spacetime.

    Parameters:
    r     : radial coordinate
    theta : polar angle (in radians)
    M     : mass of the black hole
    B     : magnetic field strength (or electric analogue)
    E     : conserved energy per unit mass
    L     : conserved angular momentum per unit mass

    Returns:
    value of V_{E,L}
    """
    sin_theta = np.sin(theta)
    Lambda = 1 + (1/4) * B**2 * r**2 * sin_theta**2
    Lambda2 = Lambda**2
    Lambda4 = Lambda2**2

    prefactor = r**4 / (L**2 * Lambda4)
    term1 = E**2 / Lambda4
    term2 = (1 - 2*M/r) * (L**2 / r**2 + 1 / Lambda2)

    return - prefactor * (term1 - term2) / 2


def Double_derivative(r, M, B, L, E):
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
        )
    )
    denominator = L**2 * (B**2 * r**2 + 4)**10
    return numerator / denominator

def E_init(r, theta, B, M):
    Lmbd = Lambda(r, theta, B)
    num = 2 * r - (r**3 * B**2 / Lmbd) - (2 * M * r**2 * B**2 / Lmbd) - (4 * M / Lmbd**2)
    den = (2 * r / Lmbd**2) - (2 * r**3 * B**2 / Lmbd**3) - (2 * M / Lmbd**2)/(1 - 2 * M / r)
    try:
        val = num / den
        if val <= 0:
            return None
        return math.sqrt(val)
    except:
        return None

def L_init(r, theta, B, M, E):
    if E is None:
        return None
    Lmbd = Lambda(r, theta, B)
    term1 = r / Lmbd
    term2 = E**2 / (Lmbd**2 * (1 - 2 * M / r)) - 1
    if term2 <= 0:
        return None
    return term1 * math.sqrt(term2)

# plt.figure()

for B in [0.01,0.05, 0.10]:
    print(f"\nFor B = {B}M:")
    r_3 = calculate_r_3(B, M, theta)
    r_1 = calculate_r_1(B, M, theta)
    print(f"Inner circular orbit: r_3 = {r_3:.4f}")
    print(f"Outer circular orbit: r_1 = {r_1:.4f}")

    V2 = []
    r_valid = []

    for r in r_range:
        E = E_init(r, theta, B, M)
        if E is None:
            continue
        L = L_init(r, theta, B, M, E)
        if L is None:
            continue
        V_val2 = Double_derivative(r, M, B, L, E)
        
        
        r_valid.append(r)
        V2.append(V_val2)

    plt.plot(r_valid, V2, label=(
        rf"$V''$ for $B = {B:.2f}M$, "
        rf"$r_3 = {r_3:.2f}M$."
    ))
    plt.axvline(r_1, color='red', linestyle='--', label=rf"$r_1 = {r_1:.2f}M$, for $B = {B}/M$.")
plt.axhline(0, color='black', linestyle='--')
plt.xlabel("r (M)",fontsize=12)
plt.ylabel(r"$V_{E,L}''$",fontsize=12)
# plt.title("Second derivative of effective potential for various B",fontsize=14)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.figure()
plt.show()


# B = 0.1  # Define B explicitly here
for B in B_list:
    r_3 = calculate_r_3(B, M, theta)
    r_1 = calculate_r_1(B, M, theta)
    r_range = np.linspace(r_3 + pert, r_1 - pert, 200)

    V = []
    r_valid = []
    E = 1
    L = 4.0  

    # Plotting
    # plt.figure()
    r_vals = np.linspace(2.5, 50, 1000)
    V_vals = potential(r_vals,theta,M,B,E,L)
    # plt.figure(figsize=(8,5))
    plt.plot(r_vals, V_vals, label=f'$B={B}/M$')
plt.axhline(0, color='gray', linestyle='--', linewidth=0.8)
plt.xlabel('$r$',fontsize=12)
plt.ylabel('$V_{E,L}$',fontsize=12)
# plt.title('Effective Potential for Timelike Observer (away from circular orbit)',fontsize=14)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.axhline(0, color='black', linestyle='--')
plt.ylim(-300,1100)
plt.show()