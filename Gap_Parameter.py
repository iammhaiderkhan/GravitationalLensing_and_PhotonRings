import math
import numpy as np
import matplotlib.pyplot as plt
from Functions import Lambda, calculate_r_1, calculate_r_3

M=1
B=0.01*M
theta=math.pi/2

# print(f"For B = {B}M")
if B:
    r_3=calculate_r_3(B,M,theta)
    r_1=calculate_r_1(B,M,theta)
    # print(f"Inner circular orbit: r_3 = {r_3:.4f}")
    # print(f"Outer circular orbit: r_1 = {r_1:.4f}")

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


def E_init(r,theta,B,M):
    Lmbd=Lambda(r,theta,B)
    num = 2 * r - (r**3 * B**2 / Lmbd) - (2 * M * r**2 * B**2 / Lmbd) - (4 * M / Lmbd**2)
    den = (2 * r / Lmbd**2) - (2 * r**3 * B**2 / Lmbd**3) - (2 * M / Lmbd**2)/(1-2*M/r)
    # print(f"For r ={r:.4f}, Energy = {num/den:.4f}")
    return math.sqrt(num/den)


def L_init(r,theta,B,M,E):
    Lmbd=Lambda(r,theta,B)
    term1 = r / Lmbd
    term2= E**2/(Lmbd**2 * (1-2*M/r)) - 1
    # print(f"For r ={r:.4f}, Energy = {E:.4f} and term2 = {term2:.4f}")
    return term1 * math.sqrt(term2)

def E_init_Schwarzschild(r, M):
    if r <= 3 * M:
        raise ValueError("No stable circular timelike orbit possible for r ≤ 3M")
    return math.sqrt((1 - 2*M/r) / (1 - 3*M/r))

def L_init_Schwarzschild(r, M):
    if r <= 3 * M:
        raise ValueError("No stable circular timelike orbit possible for r ≤ 3M")
    return math.sqrt(M * r / (1 - 3*M/r))


def eta(r,M):
    return 1-3*M/r

def b_n(n,r_s,r_O,M):
    """Compute the impact parameter b_n for given n, M, r_s (source), r_O (observer)."""
    # First, define common terms
    sqrt3 = np.sqrt(3)
    exp_term = np.exp(-(n + 0.5) * np.pi)
    
    # Denominator terms
    denom_s = 2 + (3 * M / r_s) + np.sqrt(3 + 18 * M / r_s)
    denom_O = 2 + (3 * M / r_O) + np.sqrt(3 + 18 * M / r_O)
    
    # Main expression
    prefactor = 3 * sqrt3 * M
    correction = 216 * (1 - (3 * M / r_s)) * (1 - (3 * M / r_O)) * exp_term / (denom_s * denom_O)
    
    return prefactor * (1 + correction)

def b_n_Schwarzschild(M, r_s, n):
    # Constants
    sqrt3 = np.sqrt(3) * M
    
    # Terms inside the bracket
    term1 = 1 - (3 * M / r_s)
    term2 = 2 + (3 * M / r_s) + np.sqrt(3 + 18 * M / r_s)
    term3 = 216 / (2 + sqrt3) * np.exp(-(n + 0.5) * np.pi)
    
    # Combine everything
    b_n = 3 * sqrt3 * M * (1 + (term1 / term2) * term3)
    
    return b_n

def gap_parameter(n,r_s,r_O,M):
    num = b_n(n,r_s,r_O,M)-b_n(n+1,r_s,r_O,M)
    den = b_n(n,r_s,r_O,M)
    return num/den

V2=[]
pert=0.00001
if B:
    r_range=np.linspace(r_3+pert,r_1-pert,2000)
else:
    r_range=np.linspace(3+pert,50-pert,2000)
for r in r_range:
    E=E_init(r,theta,B,M)
    V2.append(Double_derivative(r,M,B,L_init(r,theta,B,M,E),E))


# plt.plot(r_range,V)
# plt.plot(r_range,np.zeros(len(r_range)))

ISCO,OSCO = 0,0
for i in range(len(V2)):
    if V2[i-1]>0 and V2[i]<0:
        OSCO=r_range[i]
    if V2[i-1]<0 and V2[i]>0:
        ISCO=r_range[i]
# plt.xlabel("r")
# plt.ylabel(r"$V_{E,L}''$")
print(f"ISCO = {ISCO:.4f}M.")
# print(f"OSCO = {OSCO:.4f}M.") 

r_O=5000*M
if B:
    for r_s in [ISCO+pert,r_1-pert]:
        print(f"For r_s = {r_s:.2f}M and r_O = {r_O}M, b_2={b_n(2,r_s,r_O,M):.6f}, b_3={b_n(3,r_s,r_O,M):.6f} and the Gap Parameter (n={2}) is {gap_parameter(2,r_s,r_O,M):.6f}")
else: 
    for r_s in np.linspace(ISCO+pert,1500-pert,7):
        print(f"For r_s = {r_s:.2f}M and r_O = {r_O}M, b_2={b_n(2,r_s,r_O,M):.6f}, b_3={b_n(3,r_s,r_O,M):.6f} and the Gap Parameter (n={2}) is {gap_parameter(2,100,r_O,M):.6f}")


B_vals = np.linspace(0.01, 0.189366386, 20)
B_vals = np.insert(B_vals,0,0)  # Add Schwarzschild case at the start
gap_lines = []

plt.figure(figsize=(10, 6))
for idx, B in enumerate(B_vals):
    if B==0:
        ISCO = 6 * M
        # Evaluate gap parameter for n = 2, 3
        b2_ISCO = b_n_Schwarzschild(M, ISCO, 2)
        b3_ISCO = b_n_Schwarzschild(M, ISCO, 3)
        gp_ISCO = (b2_ISCO - b3_ISCO) / b2_ISCO

        b2_rO = b_n_Schwarzschild(M, 100000000, 2)
        b3_rO = b_n_Schwarzschild(M, 100000000, 3)
        gp_rO = (b2_rO - b3_rO) / b2_rO

        y_level = idx
        plt.plot([gp_ISCO, gp_rO], [y_level, y_level], label=f"B={B:.3f}", linewidth=2)
        plt.scatter([gp_ISCO, gp_rO], [y_level, y_level], color='black', s=10)
        print(f"B={B:.2f}, ISCO: {ISCO:.4f}, GP_ISCO = {gp_ISCO:.6f}, GP_rO = {gp_rO:.6f}")
    else:
        r_1 = calculate_r_1(B, M, theta)
        r_3 = calculate_r_3(B, M, theta)

        # Find ISCO from second derivative
        V = []
        pert = 1e-5
        r_range = np.linspace(r_3 + pert, r_1 - pert, 2000)
        for r in r_range:
            E = E_init(r, theta, B, M)
            L = L_init(r, theta, B, M, E)
            V.append(Double_derivative(r, M, B, L, E))

        ISCO = None
        for i in range(1, len(V)):
            if V[i - 1] < 0 and V[i] > 0:
                ISCO = r_range[i]
                break  # take the first crossing

        if ISCO is None:
            continue  # skip this B if no ISCO found

        # Evaluate gap parameter at ISCO and r_1
        r_O = 5000 * M
        gp_ISCO = gap_parameter(2, ISCO, r_O, M)
        gp_r1 = gap_parameter(2, r_1, r_O, M)

        y_level = idx  # vertical level just for visualization

        plt.plot([gp_ISCO, gp_r1], [y_level, y_level], label=f"B={B:.3f}", linewidth=2)
        # print(f"B={B:.2f}, ISCO: {ISCO:.4f}, GP = {gp_ISCO:.4f}")
        plt.scatter([gp_ISCO, gp_r1], [y_level, y_level], color='black', s=10)

plt.xlabel(r"$\Delta_2$",fontsize=12)
plt.ylabel(r"Magentic Field Strength $B$",fontsize=12)
# plt.title("Gap parameter lines between ISCO and $r_1$ for different B",fontsize=14)
plt.yticks(ticks=range(len(B_vals)), labels=[f"{B:.3f}" for B in B_vals])
# plt.axvline(gap_parameter(2,6,r_O,M), color="#3e314b", linestyle="--", linewidth=1.5)
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

for i in np.linspace(0.15,0.17,100):
    if calculate_r_1(i,M,math.pi/2)<6:
        print(f"For B = {i:.3f}/M, the overlap is none with Schwarzschild.")
        break