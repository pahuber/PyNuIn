import numpy as np
import sympy as sym
from math import factorial

        
# method to calculate radial part of Zernike polynomial
def r(n, m, rho):
    
    # if n-m even
    if (n-m) % 2 == 0:
        result = 0
        for k in np.arange(0, (n-m)/2+1, 1):
            result += ((-1)**k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*rho**(n-2*k)
        return result
    # if n-m odd
    elif (n-m) % 2 != 0:
        return 0
    # if rho is one
    elif rho == 1:
        return 1


# method to calculate full Zernike polynomial, i. e. radial plus angular part
def z(n, m, rho=sym.Symbol("rho"), theta=sym.Symbol("theta")):
    
    # check if valid input
    if n < 0 or m > n or abs(m) > n:
        print("Make sure that n >= 0 and |m| <= n")
    else:
        if m >= 0:
            return r(n, m, rho)*sym.cos(m*theta)
        else:
            return r(n, abs(m), rho)*sym.sin(abs(m)*theta)


# method to simulate wavefront error (WFE) by summing several Zernike polynomials
def wfe_basic(rho=sym.Symbol("rho"),
        theta=sym.Symbol("theta"),
        piston=0,
        tilty=0,
        tiltx=0,
        astig_oblique=0,
        defocus=0,
        astig_vertical=0):
    
    return piston * z(0, 0, rho, theta) + tilty * z(1, -1, rho, theta) + tiltx * z(1, 1, rho, theta) + astig_oblique * z(2, -2, rho, theta) + defocus * z(2, 0, rho, theta) + astig_vertical * z(2, 2, rho, theta)


# method to simulate general wavefront error
def wfe(indices_coeffs,
        rho=sym.Symbol("rho"),
        theta=sym.Symbol("theta"),):
    
    summation = 0
    for el in indices_coeffs:
        n = el[0]
        m = el[1]
        coeff = el[2]
        summation += coeff * z(n, m, rho, theta)
        
    return summation


# method to simulate an ideal wavefront amplitude
def amp_id(a0, t=sym.Symbol("t"), lam=450e-10, c=3e8):
    return a0 * np.exp(2j*np.pi*c*t/lam)


# method to simulate an aberrated wafefront amplitude, requires the definition of a wavefront error
def amp_ab(a0, wfe, t=sym.Symbol("t"), lam=450e-10, rho=sym.Symbol("rho"), theta=sym.Symbol("theta"), c=3e8):
    return a0 * np.exp(2j*np.pi/lam*(c*t-wfe))




