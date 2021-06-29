import numpy as np
import sympy as sym
from math import factorial
from pynuin.util.general import kronecker_delta
from scipy.optimize import fsolve
from scipy.integrate import dblquad


def r(n, m, r, r_max):
    '''
    Returns the radial part of Zernike polynomials.

            Parameters:
                    n (int): First index of Zernike polynomial
                    m (int): Second index of Zernike polynomial
                    r (float): Radial coordinate
                    r_max (float): Maximum of radial coordinate, such that rho=r/r_max in [0, 1]

            Returns:
                    (float): Radial part of Zernike polynomial
    '''

    # if n-m even
    if (n-m) % 2 == 0:

        # if rho is one
        # if r == r_max:
        #     return 1

        # else:
        result = 0
        for k in np.arange(0, (n-m)/2+1, 1):
            result += ((-1)**k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*(r/r_max)**(n-2*k)
        return result
    # if n-m odd
    elif (n-m) % 2 != 0:
        return 0


# method to calculate full Zernike polynomial, i. e. radial plus angular part
def z(n, m, rho, theta, rho_max=1):
    '''
    Returns the full Zernike polynomial Z_n^m.

            Parameters:
                    n (int): First index of Zernike polynomial
                    m (int): Second index of Zernike polynomial
                    r (float): Radial coordinate
                    theta (float): Angular coordinate
                    r_max (float): Maximum of radial coordinate, such that rho=r/r_max in [0, 1]

            Returns:
                    (float): Zernike polynomial
    '''

    # check if valid input
    if n < 0 or m > n or abs(m) > n:
        print("Make sure that n >= 0 and |m| <= n")
    else:
        sign = np.sign(m)
        if sign == 0:
            sign = 1
        norm = sign*np.sqrt((2*(n+1))/(1+kronecker_delta(m, 0)))
        if m >= 0:
            return norm*r(n, m, rho, rho_max)*np.cos(m*theta)
        else:
            return norm*r(n, abs(m), rho, rho_max)*np.sin(abs(m)*theta)


# method to simulate general wavefront error
def wfe(indices_coeffs,
        rho,
        theta,
        rho_max):
    '''
    Returns a wavefront error composed of a sum of several Zernike polynomial terms Z_n^m.

            Parameters:
                    indices_coeffs (list): List containing tuples of the kind (n, m, coeff)
                    r (float): Radial coordinate
                    theta (float): Angular coordinate
                    r_max (float): Maximum of radial coordinate, such that rho=r/r_max in [0, 1]

            Returns:
                    (float): Wavefront error
    '''

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
    '''
    Returns the coefficient z_n^m = z of a sum of Zernike polynomials assuming the same coefficient for each contribution.

            Parameters:
                    rms (float): Total root mean square wavefront aberration
                    indices (list): List containing tuples of the kind (n, m)

            Returns:
                    (float): Zernike polynomial coefficient
    '''

    return abs(rms/np.sqrt(len(indices)))


def noll_index(n, m):
    j = (n*(n+1))/2+abs(m)

    if m > 0 and (n % 4 == 0 or n % 4 == 1):
        j += 0
    elif m < 0 and (n % 4 == 2 or n % 4 == 3):
        j += 0
    elif m >= 0 and (n % 4 == 2 or n % 4 == 3):
        j += 1
    elif m <= 0 and (n % 4 == 0 or n % 4 == 1):
        j += 1

    return j
