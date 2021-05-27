import numpy as np
import sympy as sym
from math import factorial
        

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
        if m >= 0:
            return r(n, m, rho, rho_max)*sym.cos(m*theta)
        else:
            return r(n, abs(m), rho, rho_max)*sym.sin(abs(m)*theta)
        

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