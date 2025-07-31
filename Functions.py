#Includes most of the functions that are repeatedly used in different files

import numpy as np
import math
from scipy.integrate import quad

# Define the Lambda function (depends on r and theta)
def Lambda(r, theta, B):
    return 1 + (1/4) * B**2 * r**2 * np.sin(theta)**2

# Define the system of ODEs
def equations(lambda_param, y, M, B, pt, pphi):
    """
    y = [t, phi, r, theta, pr, ptheta]
    lambda_param: Affine parameter
    M: Mass of the black hole
    B: Magnetic field strength
    pt: Conserved energy
    pphi: Conserved angular momentum
    """
    t, phi, r, theta, pr, ptheta = y

    # Lambda value
    L = Lambda(r, theta, B)

    # Equations for the coordinates
    dt_dlambda = -pt / (L**2 * (1 - 2 * M / r))
    dphi_dlambda = pphi * L**2 / (r**2 * np.sin(theta)**2)
    dr_dlambda = pr * (1 - 2 * M / r) / L**2
    dtheta_dlambda = ptheta / (r**2 * L**2)

    # Equations for the momenta
    dpr_dlambda = (
        (2 * L**2 * M / (r - 2 * M)**2 - L * r**2 * B**2 * np.sin(theta)**2) * dr_dlambda**2
        - (2 * L**2 * r + L * r**3 * B**2 * np.sin(theta)**2) * dtheta_dlambda**2
    )

    dptheta_dlambda = -(
        (1 / (1 - 2 * M / r) * dr_dlambda**2 + r**2 * dtheta_dlambda**2) * L * r**2 * B**2 * np.sin(theta) * np.cos(theta)
    )

    return [dt_dlambda, dphi_dlambda, dr_dlambda, dtheta_dlambda, dpr_dlambda, dptheta_dlambda]

#Defining initial value of pt in accordance with the null condition for circular orbit
#Need to redefine this function for not a circular orbit
def pt_init(pphi,theta,M,B,r):
    return -pphi*Lambda(r,theta,B)**2/(r*math.sin(theta))*(1-2*M/r)**(1/2)

def L_init(r,theta,B,M,E):
    # Compute numerator and denominator of the L^2/E^2 equation
    numerator = r**3 * np.sin(theta)**2 * (r**2 * B**2 * np.sin(theta)**2 * (r - 2*M) + 2*M * Lambda(r, theta, B))
    denominator = Lambda(r, theta, B)**4 * (r - 2*M)**2 * (2*Lambda(r, theta, B) - r**2 * B**2 * np.sin(theta)**2)

    # Compute L
    if denominator <= 0:
        raise ValueError("Invalid parameters leading to non-physical L^2/E^2.")

    L_squared = (E**2) * (numerator / denominator)
    L = np.sqrt(L_squared) if L_squared > 0 else np.nan  # Ensure non-negative L

    return L
    
def pphi_init(pt, pr, ptheta, theta, M, B, r):
    lambd=Lambda(r,theta,B)
    pphi=r*math.sin(theta)/lambd*(pt**2/(lambd**2*(1-2*M/r))-pr**2*(1-2*M/r)/lambd**2-ptheta**2/(r**2*lambd**2))**0.5
    return pphi


def calculate_r_1(B, M, theta):
    # Compute the discriminant term under the square root
    discriminant = -375 * (B*math.sin(theta))**10 * M**4 + 1352 * (B*math.sin(theta))**8 * M**2 - 48 * (B*math.sin(theta))**6

    # Handle the square root of the discriminant
    if discriminant < 0:
        sqrt_term = np.sqrt(3) * np.sqrt(abs(discriminant)) * 1j  # Handle complex case
    else:
        sqrt_term = np.sqrt(3) * np.sqrt(discriminant)

    # Compute the main term argument
     # Compute the main term argument
    term1_argument = 125 * (B*math.sin(theta))**6 * M**3 - 1188 * (B*math.sin(theta))**4 * M + 18 * sqrt_term

    # Function to compute the cube root correctly (handles real and complex numbers)
    def cube_root(x):
        return np.sign(x) * abs(x)**(1/3) if np.isreal(x) else x**(1/3)

    # Compute the cube root of the main term
    term1_cbrt = cube_root(term1_argument)
    term2 = (-25 * (B*math.sin(theta))**4 * M**2 - 36 * (B*math.sin(theta))**2) / (9 * (B*math.sin(theta))**2 * term1_cbrt)

    r = term1_cbrt / (9 * (B*math.sin(theta))**2) - term2 + (5 * M) / 9

    # Discard tiny imaginary parts (if any)
    if abs(r.imag) < 1e-10:
        r = r.real  # Keep only the real part
    return r

def calculate_r_2(B, M, theta):
    """
    Calculate the value of r using the given analytical expression.

    Parameters:
        B (float): Magnetic field strength.
        M (float): Black hole mass (default is 1.0).

    Returns:
        r (float): The computed value of r (with tiny imaginary parts discarded).
    """
    # Compute the discriminant term under the square root
    # Compute the discriminant term under the square root
    discriminant = -375 * (B*math.sin(theta))**10 * M**4 + 1352 * (B*math.sin(theta))**8 * M**2 - 48 * (B*math.sin(theta))**6

    # Handle the square root of the discriminant
    if discriminant < 0:
        sqrt_term = np.sqrt(3) * np.sqrt(abs(discriminant)) * 1j  # Handle complex case
    else:
        sqrt_term = np.sqrt(3) * np.sqrt(discriminant)

    # Compute the main term argument
    term1_argument = 125 * (B*math.sin(theta))**6 * M**3 - 1188 * (B*math.sin(theta))**4 * M + 18 * sqrt_term

    # Function to compute the cube root correctly (handles real and complex numbers)
    def cube_root(x):
        return np.sign(x) * abs(x)**(1/3) if np.isreal(x) else x**(1/3)

    # Compute the cube root of the main term
    term1_cbrt = cube_root(term1_argument)

    # Compute the two terms in the equation
    term1 = -((1 - 1j * np.sqrt(3)) * term1_cbrt) / (18 * (B*math.sin(theta))**2)
    term2 = ((1 + 1j * np.sqrt(3)) * (-25 * (B*math.sin(theta))**4 * M**2 - 36 * (B*math.sin(theta))**2)) / (18 * (B*math.sin(theta))**2 * term1_cbrt)

    # Compute the final r value
    r = term1 + term2 + (5 * M) / 9

    # Discard tiny imaginary parts (if any)
    if abs(r.imag) < 1e-10:
        r = r.real  # Keep only the real part

    return r

def calculate_r_3(B, M, theta):
    """
    Calculate the value of r using the given analytical expression.

    Parameters:
        B (float): Magnetic field strength.
        M (float): Black hole mass (default is 1.0).

    Returns:
        r (float): The computed value of r (with tiny imaginary parts discarded).
    """
    # Compute the discriminant term under the square root
    discriminant = -375 * (B*math.sin(theta))**10 * M**4 + 1352 * (B*math.sin(theta))**8 * M**2 - 48 * (B*math.sin(theta))**6

    # Handle the square root of the discriminant
    if discriminant < 0:
        sqrt_term = np.sqrt(3) * np.sqrt(abs(discriminant)) * 1j  # Handle complex case
    else:
        sqrt_term = np.sqrt(3) * np.sqrt(discriminant)

    # Compute the main term argument
    term1_argument = 125 * (B*math.sin(theta))**6 * M**3 - 1188 * (B*math.sin(theta))**4 * M + 18 * sqrt_term

    # Function to compute the cube root correctly (handles real and complex numbers)
    def cube_root(x):
        return np.sign(x) * abs(x)**(1/3) if np.isreal(x) else x**(1/3)

    # Compute the cube root of the main term
    term1_cbrt = cube_root(term1_argument)

    # Compute the two terms in the equation
    term1 = -((1 + 1j * np.sqrt(3)) * term1_cbrt) / (18 * (B*math.sin(theta))**2)
    term2 = ((1 - 1j * np.sqrt(3)) * (-25 * (B*math.sin(theta))**4 * M**2 - 36 * (B*math.sin(theta))**2)) / (18 * (B*math.sin(theta))**2 * term1_cbrt)

    # Compute the final r value
    r = term1 + term2 + (5 * M) / 9

    # Discard tiny imaginary parts (if any)
    if abs(r.imag) < 1e-10:
        r = r.real  # Keep only the real part

    return r


def V_eff(M, B, E, L, r, theta):
    """
    Computes the effective potential V_eff for given values of M, B, E, L, and r.
    """
    Lambda = (1 + (r**2 * (B*math.sin(theta))**2) / 4)
    return (-E**2/(Lambda**2*(1-2*M/r))+L**2*Lambda**2/(r**2*(math.sin(theta))**2))
    
def V_eff_first_derivative(M, B, E, L, r, theta):
    # Define Lambda function
    Lambda = 1 + (r**2 * (B*np.sin(theta))**2) / 4

    # Compute the terms
    term1 = E**2 * ((r**2 * B**2 * np.sin(theta)**2) / (Lambda**3 * (r - 2 * M)) + (2 * M) / (Lambda**2 * (r - 2 * M)**2))
    term2 = L**2 * ((Lambda * r * B**2 * np.sin(theta)**2) / (r**2 * np.sin(theta)**2) - (2 * Lambda**2) / (r**3 * np.sin(theta)**2))

    # Sum the terms
    result = term1 + term2

    return result

def V_eff_second_derivative(M, B, E, L, r, theta):
    """
    Evaluates the given second derivative expression for given values of M, B, E, L, and r.
    """
    # Define common terms
    Lambda = (1 + (r**2 * (B*math.sin(theta))**2) / 4)
    denom1 = Lambda**2 * (1 - (2 * M) / r)

    # Compute the right-hand side
    term1 = - (4 * (B*math.sin(theta))**2 * L**2 * Lambda) / r**2
    term2 = (6 * L**2 * Lambda**2) / r**4
    term3 = (L**2 * (((B*math.sin(theta))**4 * r**2) / 2 + (B*math.sin(theta))**2 * Lambda)) / r**2

    term4 = - E**2 * (((8 * M**2) / (r**4 * (1 - (2 * M) / r)**3) + (4 * M) / (r**3 * (1 - (2 * M) / r)**2))) / Lambda**2
    term5 = - E**2 * (4 * (B*math.sin(theta))**2 * M) / (r * Lambda**3 * (1 - (2 * M) / r)**2)
    term6 = - E**2 * ((3 * (B*math.sin(theta))**4 * r**2) / (2 * Lambda**4) - (B*math.sin(theta))**2 / Lambda**3) / (1 - (2 * M) / r)

    rhs = term1 + term2 + term3 + term4 + term5 + term6

    return rhs

def V_eff_first_derivative_theta(M, B, E, L, r, theta):
    """
    Calculate the value of first derivative of V_eff wrt theta

    Parameters:
        B (float): Magnetic field strength.
        M (float): Black hole mass (default is 1.0).
        E (float): Conserved energy
        L (float): Conserved angular momentum (z-component)
        r (float): the radial value at which the V_eff is being determined
        theta (float): the theta value at which the V_eff is being determined

    Returns:
        
    """
    term1=(E**2*math.sin(theta)*math.cos(theta)*r**3*B**2)/((Lambda(r,theta,B)**3)*(r-2*M))
    term2=(2*L**2*(Lambda(r,theta,B)**2)*math.cos(theta))/(r**2*(math.sin(theta))**3)
    term3=(L**2*Lambda(r,theta,B)*B**2*math.cos(theta))/(math.sin(theta))
    return term1-term2+term3

def V_eff_first_derivative_theta2(M, B, E, L, r, theta):
    """
    Calculate the value of first derivative of V_eff wrt theta (without the cos theta factor)

    Parameters:
        B (float): Magnetic field strength.
        M (float): Black hole mass (default is 1.0).
        E (float): Conserved energy
        L (float): Conserved angular momentum (z-component)
        r (float): the radial value at which the V_eff is being determined
        theta (float): the theta value at which the V_eff is being determined

    Returns:
        
    """
    term1=(E**2*math.sin(theta)*r**3*B**2)/((Lambda(r,theta,B)**3)*(r-2*M))
    term2=(2*L**2*(Lambda(r,theta,B)**2))/(r**2*(math.sin(theta))**3)
    term3=(L**2*Lambda(r,theta,B)*B**2)/(math.sin(theta))
    return term1-term2+term3

import numpy as np

def V_eff_second_derivative_theta(M, B, E, L, r, theta):
    # Compute individual terms
    term1 = L**2 * (2/r**2 - (B**4 * r**2) / 8)
    term2 = -(64 * E**2 * B**2 * r**3) / ((B**2 * r**2 + 4)**3 * (r - 2*M))
    # Combine terms
    result = term1 + term2
    return result


def p(B,M, theta):
    return (-36-25*M*B**2*(math.sin(theta))**2)/(27*B**2*(math.sin(theta))**2)

def q(B,M,theta):
    return (-250*M**3*B**2*(math.sin(theta))**2+2376*M)/(729*B**2*(math.sin(theta))**2)

def Disc(p,q):
    return -4*p**3-27*q**2

def cubic_equation(r,B,M,theta):
    return (3*(B*math.sin(theta))**2)*r**3-5*(M*(B*math.sin(theta))**2)*r**2-4*r+12*M

# For the bounded light-like regions
def dr_dlambda(lambda_param,r,E,B,L,M,sign):
    """
    Calculate the value of first derivative of r wrt the affine parameter lambda

    Parameters:
        B (float): Magnetic field strength.
        M (float): Black hole mass (default is 1.0).
        E (float): Conserved energy
        L (float): Conserved angular momentum (z-component)
        r (float): the radial value at which the V_eff is being determined
        sign (int): +1 or -1 indicating the direction of the lightray (negative for ingoing)

    Returns:
        The instantaneous value of the derivative
    """
    term1=E**2/(1+r**2*B**2/4)**4
    term2=-L**2/(r**2*(1-2*M/r))
    return sign*np.sqrt(term1-term2)

def dphi_dlambda(lambda_param,r,B,L):
    """
    Calculate the value of first derivative of phi wrt the affine parameter lambda

    Parameters:
        B (float): Magnetic field strength.
        L (float): Conserved angular momentum (z-component)
        r (float): the radial value at which the V_eff is being determined

    Returns:
        The instantaneous value of the derivative
    """
    return L*(1+r**2*B**2/4)**2/r**2

def H(r,B, M, theta):
    return Lambda(r,theta,B)**2 * math.sqrt(1-2*M/r)/(r * math.sin(theta))

def dH_dr(r, B, M, theta):
    term1=Lambda(r,theta,B) * B**2 * math.sin(theta) * math.sqrt(1-2*M/r)
    term2=-Lambda(r,theta,B)**2 * math.sqrt(1-2*M/r) / (r**2 * math.sin(theta))
    term3=Lambda(r,theta,B)**2 * M / (r**3 * math.sin(theta) * math.sqrt(1-2*M/r))

    return term1+term2+term3

def dH_dtheta(r, B, M, theta):
    term1=Lambda(r,theta,B) * r * B**2 * math.cos(theta) * math.sqrt(1-2*M/r)
    term2=-Lambda(r,theta,B)**2 * math.sqrt(1-2*M/r) *math.cos(theta) / (r * math.sin(theta)**2)
    return term1+term2

import numpy as np

def theta_dot_squared(r, theta, B, M, E, L):
    # Compute Lambda
    Lambda = 1 + 0.25 * B**2 * r**2 * np.sin(theta)**2

    # Compute terms step by step
    denom = Lambda * r**3 * B**2 * np.sin(theta)**2 + 2 * r * Lambda**2
    term1 = (r * B**2 * np.sin(theta)**2) / (Lambda**3 * (1 - 2 * M / r))
    term2 = (2 * M) / (r**2 * Lambda**2 * (1 - 2 * M / r)**2)
    term3 = (B**2 * Lambda / r) - (2 * Lambda**2) / (r**3 * np.sin(theta)**2)

    # Final expression
    result = (1 / denom) * ((term1 + term2) * E**2 + term3 * L**2)

    return result

import numpy as np

def theta_dot_squared2(r, theta, M, B, L):
    """
    Evaluates the expression for \dot{\theta}^2 for lightlike geodesics in the Ernst spacetime.

    Parameters:
    r     : radial coordinate
    theta : polar angle (in radians)
    M     : mass parameter
    B     : magnetic field parameter
    L     : angular momentum

    Returns:
    dot_theta_sq : the value of \dot{\theta}^2
    """
    sin_theta = np.sin(theta)
    sin2 = sin_theta**2
    Lambda = 1 + 0.25 * B**2 * r**2 * sin2
    one_minus_2M_over_r = 1 - 2*M/r

    # Precompute repetitive terms
    A = r * B**2 * sin2 / (Lambda**3 * one_minus_2M_over_r)
    B_term = 2*M / (r**2 * Lambda**2 * one_minus_2M_over_r**2)
    C = Lambda**4 * one_minus_2M_over_r / (r**2 * sin2)

    numerator = (A + B_term) * C + (B**2 * Lambda / r) - (2 * Lambda**2 / (r**3 * sin2))

    D = Lambda * r**3 * B**2 * sin2
    E = 2 * r * Lambda**2
    F = (A + B_term) * Lambda**4 * one_minus_2M_over_r * r**2

    denominator = D + E - F

    dot_theta_sq = (numerator / denominator) * L**2
    return dot_theta_sq

def impact_parameter(M,r_m,B):
    return r_m / ((1 + 1/4 * r_m**2 * B**2)**2 * (math.sqrt(1- 2*M/r_m)))

def integrand(r, r_m, M, B):
    # Compute the expression under the square root
    numerator = (1 - 2 * M / r_m) * r**4
    term2 = r**2 * r_m**2
    term3 = 2 * M * r * r_m**2
    denominator = np.sqrt(numerator - term2 + term3)
    
    if denominator <= 0:
        print("Negative denominator)")
        return 0  # avoid complex numbers and division by zero
    return r_m / denominator

def deflection_angle(r_m, M, B):
    integral, error = quad(integrand, r_m, np.inf, args=(r_m, M, B), limit=100)
    return 2 * integral - np.pi

def geodesic_equations(lambda_param, state, M, B, E, Q):
    r, v_r, theta, v_theta = state
    ddot_r = -((1 - 2 * M / r) / Lambda(r, theta, B)**2) * (
    (Lambda(r, theta, B) / (1 - 2 * M / r)) *
    (r * B**2 * np.sin(theta)**2 - Lambda(r, theta, B) * M / (r**2 * (1 - 2 * M / r))) * v_r**2
    + (Q * B**2 * np.sin(theta) * np.cos(theta)) / (Lambda(r, theta, B) * (1 - 2 * M / r))* v_r
    - Q / (r**3 * Lambda(r, theta, B)**2)
    + M * E**2 / (r**2 * Lambda(r, theta, B)**2 * (1 - 2 * M / r)**2)
)


    ddot_theta = -2 * Q * (
        (v_r / (r**3 * Lambda(r, theta, B)**2)) +
        (1 / (r**2 * Lambda(r, theta, B)**3)) *
        (r * B**2 * np.sin(theta)**2 * v_r +
        r**2 * B**2 * np.sin(theta) * np.cos(theta) * v_theta)
    )
    return [v_r, ddot_r, v_theta, ddot_theta]
