import numpy as np
import sympy as sym
from math import factorial
from pynuin.util.general import kronecker_delta
from scipy.optimize import fsolve
from scipy.integrate import dblquad
        

# method to calculate radial part of Zernike polynomial
def r(n, m, rho, rho_max):
    
    # if n-m even
    if (n-m) % 2 == 0:
        result = 0
        for k in np.arange(0, (n-m)/2+1, 1):
            result += ((-1)**k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*(rho/rho_max)**(n-2*k)
        return result
    # if n-m odd
    elif (n-m) % 2 != 0:
        return 0
    # if rho is one
    elif rho == 1:
        return 1
    

# method to calculate full Zernike polynomial, i. e. radial plus angular part
def z(n, m, rho, theta, rho_max=1):
    
    # check if valid input
    if n < 0 or m > n or abs(m) > n:
        print("Make sure that n >= 0 and |m| <= n")
    else:
        norm = np.sign(m)*np.sqrt((2*(n+1))/(1+kronecker_delta(m, 0)))
        if m >= 0:
            return norm*r(n, m, rho, rho_max)*sym.cos(m*theta)
        else:
            return norm*r(n, abs(m), rho, rho_max)*sym.sin(abs(m)*theta)
        

# method to simulate general wavefront error
def wfe(indices_coeffs,
        rho,
        theta,
        rho_max=1):
    
    summation = 0
    for el in indices_coeffs:
        n = el[0]
        m = el[1]
        coeff = el[2]
        summation += coeff * z(n, m, rho, theta, rho_max)
        
    return summation


# def get_coeff_from_rms(rms, indices):
    
#     def integrand1(rho, theta, coeff, indices):        
#         wfe = 0
#         for el in indices:
#             n = el[0]
#             m = el[1]
#             wfe += coeff * z(n, m, rho, theta)
            
#         return wfe**2*rho
    
#     def integrand2(rho, theta, coeff, indices):        
#         wfe = 0
#         for el in indices:
#             n = el[0]
#             m = el[1]
#             wfe += coeff * z(n, m, rho, theta)
            
#         return wfe*rho
    
#     def func(coeff):
#         part1 = 1/np.pi* dblquad(integrand1, 0, 2*np.pi, lambda rho: 0, lambda rho:1, args=(coeff, indices))[0]
#         part2 = 1/np.pi**2* (dblquad(integrand2, 0, 2*np.pi, lambda rho: 0, lambda rho:1, args=(coeff, indices))[0])**2
#         function = np.sqrt(part1 - part2)
#         return function - rms
    
#     coeff = fsolve(func, 1e-5)[0]
    
#     return coeff


def get_coeff_from_rms(rms, indices):
    
    return abs(rms/np.sqrt(len(indices)))