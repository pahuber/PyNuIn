import numpy as np
import matplotlib.pyplot as plt
from zernike import RZern
import poppy as pp
from matplotlib import cm
from mpl_toolkits import mplot3d
import sympy as sym
from math import factorial
import pylab as pb





#define Zernike polynomials

        
# method to calculate radial part of Zernike polynomial
def R(n, m, rho):
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
# TODO: check arguments
def Z(n, m, rho=sym.Symbol("rho"), theta=sym.Symbol("theta")):
    # check if valid input
    if n < 0 or m > n:
        print("Make sure that n >= 0 and |m| <= n")
    else:
        if m >= 0:
            return R(n, m, rho)*sym.cos(m*theta)
        else:
            return R(n, abs(m), rho)*sym.sin(abs(m)*theta)


def WFE(rho=sym.Symbol("rho"), theta=sym.Symbol("theta"), piston=0, tilty=0, tiltx=0, astig_oblique=0, defocus=0, astig_vertical=0):
    return piston * Z(0, 0, rho, theta) + tilty * Z(1, -1, rho, theta) + tiltx * Z(1, 1, rho, theta) + astig_oblique * Z(2, -2, rho, theta) + defocus * Z(2, 0, rho, theta) + astig_vertical * Z(2, 2, rho, theta)



def A_id(a0, t=sym.Symbol("t"), lam=450e-10, c=3e8):
    return a0 * np.exp(2j*np.pi*c*t/lam)


def A_ab(a0, wfe, t=sym.Symbol("t"), lam=450e-10, rho=sym.Symbol("rho"), theta=sym.Symbol("theta"), c=3e8):
    return a0 * np.exp(2j*np.pi/lam*(c*t-wfe))




