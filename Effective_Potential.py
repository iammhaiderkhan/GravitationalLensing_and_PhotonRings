# Helps plot the equatorial lightlike potential

import numpy as np 
import math
from Functions import Lambda, calculate_r_3, calculate_r_1
import matplotlib.pyplot as plt


x_range=10000
x=np.linspace(2, 35,x_range)
theta=math.pi/2

M = 1
r3=calculate_r_3(0.12,M,theta)
r1=calculate_r_1(0.12,M,theta)


def V_eff(B,r,M):
    term1 = Lambda(r,theta,B)**4/r**2
    term2 = (1 - 2*M/r)

    return term1*term2
    

#   Plotting V_eff vs r
B_list=[0.08, 0.12, 0.189366386, 0.25, 0.3]
for B in B_list:

    # E=pt_init(L,theta,M,B,r)
    # L = 3 * np.sqrt(3) * M
    x_axis=[]
    for i in range(x_range):
        x_axis.append(0)

    V2=[]
    for i in x:
        V2.append(V_eff(B,i,M))
    if B == 0.189366386:
        plt.plot(x,V2,label=rf"$B = B_c$")
    else:
        plt.plot(x,V2,label=rf"$B = ${B}/M")
    plt.xlabel("r", fontsize=12)
    plt.ylabel(r"$V_{eff}$",fontsize=12)
    plt.ylim(0,0.15)
    # plt.xlim(2,3)
    # plt.title(f"V_eff vs r for B={B}/M")
    plt.plot(x,x_axis)

# r_s=17

# plt.axhline(V_eff(B,r3,M), color="#151444", linestyle="--", linewidth=1.5, label=r"$V_{eff}(r_2)$")
# plt.axhline(V_eff(B,r_s,M), color="#7a5e96", linestyle="--", linewidth=1.5, label=r"$\frac{E^2}{L^2}|_{crit}$")
# plt.axvline(r3, color="#020507", linestyle=":", linewidth=2, label=r"$r_2$")
# plt.axvline(15.205, color="#DA0F0F", linestyle=":", linewidth=2, label=r"$r_O'$")
# plt.axvline(r_s, color="#DA0F0F", linestyle="-", linewidth=2, label=r"$r_S>r_O'$")
    plt.legend()
plt.show()