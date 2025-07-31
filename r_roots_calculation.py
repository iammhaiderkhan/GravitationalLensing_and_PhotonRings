import math
import matplotlib.pyplot as plt
from Functions import Lambda, pt_init
from Functions import calculate_r_1, calculate_r_2, calculate_r_3
from Functions import V_eff, V_eff_second_derivative, V_eff_second_derivative_theta
from Functions import p, q, Disc


# Initial conditions and parameters
M = 1.0  # Black hole mass
pphi = 2.0  # Conserved angular momentum
p_theta = 0.0  # Angular motion in theta
p_r = 0.0  # Radial momentum
B_values = [i/100 for i in range(5,20)]
theta=math.pi/2
B_c=0.189366386/(M*math.sin(theta))
B_c2=1.889305085/(M*math.sin(theta))

#Printing values of r_1, r_2 and r_3 for a range of B
#Also printing the respective stability
for B in B_values:
    r_1=calculate_r_1(B, M,theta)
    r_2=calculate_r_2(B, M,theta)
    r_3=calculate_r_3(B, M,theta)
    pt1=pt_init(pphi,theta,M,B,r_1)
    if not isinstance(r_2, complex):
        pt2=pt_init(pphi,theta,M,B,r_2)
        pt3=pt_init(pphi,theta,M,B,r_3)
    Stability_1=""
    Stability_2=""
    Stability_3=""
    if V_eff_second_derivative_theta(M, B, pt1, pphi, r_1, theta)>0:
        Stability_1="Stable"
    else:
        Stability_1="Unstable"
    if V_eff_second_derivative_theta(M, B, pt2, pphi, r_2, theta)>0:
        Stability_2="Stable"
    else:
        Stability_2="Unstable"
    if V_eff_second_derivative_theta(M, B, pt3, pphi, r_3, theta)>0:
        Stability_3="Stable"
    else:
        Stability_3="Unstable"
    print("----------------------------------------------------------------------------------------------------------------")
    if Disc(p(B,1,theta),q(B,1,theta))>0:
        print(f"For B = {B} and theta = pi/{1/(theta/math.pi):.0f}, number of real roots: 3, r_1={r_1:.4f} ({Stability_1}), r_2={r_2:.4f}  ({Stability_2}), r_3={r_3:.4f}  ({Stability_3})")
    else:
        print(f"For B = {B} and theta = pi/{1/(theta/math.pi):.0f}, number of real roots: 1, r_1={r_1:.4f}  ({Stability_1}), r_2={r_2:.4f}  ({Stability_2}), r_3={r_3:.4f}  ({Stability_3})")
print("----------------------------------------------------------------------------------------------------------------")



r1,r2,r3,B=[],[],[],[]
for b in [i/1000 for i in range(5,2000)]:
    B.append(b)
    r1.append(calculate_r_1(b,M,theta))
    if b<B_c or b>B_c2:
        r2.append(calculate_r_2(b,M,theta))
        r3.append(calculate_r_3(b,M,theta))
    else:
        r2.append(0)
        r3.append(0)

twom=[]
x_axis=[]
for i in range(len(B)):
    x_axis.append(0)
    twom.append(2)

print(f"Critical B: {B_c:.4f}")
# for b in [i/1000 for i in range(150,200)]:
#     r=calculate_r_1(b,M,theta)
#     if V_eff_second_derivative(M,b,pt_init(pphi,theta,M,b,r),pphi,r,theta)<0:
#         print("r_1 becomes unstable for B=",b)
#         break
# for b in [i/1000 for i in range(1500,1700)]:
#     r=calculate_r_1(b,M,theta)
#     if r<2:
#         print("r_1 goes behind the horizon for B=",b)
#         break

plt.plot(B, r1, label="$r_1$")
plt.plot(B, r3, label="$r_2$")
plt.plot(B, r2, label="$r_3$")
plt.plot(B, twom, label="$r=2m$")
plt.plot(B,x_axis, label="x-axis")
plt.xlabel("B")
plt.ylabel("r")
plt.title(f"Circular light-like geodesic r values vs B for theta= pi/{1/(theta/math.pi):.0f}")
plt.ylim(-10, 10)
plt.legend()
plt.show()