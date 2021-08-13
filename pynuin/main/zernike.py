import numpy as np
from math import factorial
from pynuin.util.general import kronecker_delta
from scipy.optimize import fsolve


def n_m_from_noll(j):
    '''
    Returns two ordinary indices n, m for Zernike polynomials given a single index j according to Noll's convention.

            Parameters:
                    j (int): Noll index of Zernike polynomial

            Returns:
                    (int, int):  Ordinary indices n, m of Zernike polynomial
    '''
    
    nmax = int(j/2)
    nmin = 0
    mmin = -nmax
    mmax = nmax
    
    for n in range(nmin, nmax+1, 1):
        for m in range(mmin, mmax+1, 1):
            if noll_index(n, m) == j and (n-m) % 2 == 0 and abs(m)<=n:
                return n, m
                

def noll_index(n, m):
    '''
    Returns a single index for each pair n, m of Zernike indices according to Noll's convention.

            Parameters:
                    n (int): First index of Zernike polynomial
                    m (int): Second index of Zernike polynomial

            Returns:
                    (int): Noll index of Zernike polynomial
    '''
    
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


def r(j, r, r_max):
    '''
    Returns the radial part of Zernike polynomials.

            Parameters:
                    j (int): Noll index of Zernike polynomial
                    r (float): Radial coordinate
                    r_max (float): Maximum of radial coordinate, such that rho=r/r_max in [0, 1]

            Returns:
                    (float): Radial part of Zernike polynomial
    '''

    n, m = n_m_from_noll(j)
    m = abs(m)
    
    # if n-m even
    if (n-m) % 2 == 0:
        result = 0
        for k in np.arange(0, (n-m)/2+1, 1):
            result += ((-1)**k*factorial(n-k))/(factorial(k)*factorial((n+m)/2-k)*factorial((n-m)/2-k))*(r/r_max)**(n-2*k)
        return result
    # if n-m odd
    elif (n-m) % 2 != 0:
        return 0


def z(j, rho, theta, rho_max=1):
    '''
    Returns the full Zernike polynomial Z_n^m.

            Parameters:
                    j (int): Noll index of Zernike polynomial
                    r (float): Radial coordinate
                    theta (float): Angular coordinate
                    r_max (float): Maximum of radial coordinate, such that rho=r/r_max in [0, 1]

            Returns:
                    (float): Zernike polynomial
    '''

    n, m = n_m_from_noll(j)
    
    # check if valid input
    if n < 0 or m > n or abs(m) > n:
        print("Make sure that n >= 0 and |m| <= n")
    else:
        sign = np.sign(m)
        if sign == 0:
            sign = 1
        norm = sign*np.sqrt((2*(n+1))/(1+kronecker_delta(m, 0)))
        if m >= 0:
            return norm*r(j, rho, rho_max)*np.cos(m*theta)
        else:
            return norm*r(j, rho, rho_max)*np.sin(abs(m)*theta)


def wfe(indices_coeffs,
        rho,
        theta,
        rho_max):
    '''
    Returns a wavefront error composed of a sum of several Zernike polynomial terms Z_n^m.

            Parameters:
                    indices_coeffs (list): List containing tuples of the kind (j, coeff), where j is the Zernike poylnomial index according to Noll's convention
                    r (float): Radial coordinate
                    theta (float): Angular coordinate
                    r_max (float): Maximum of radial coordinate, such that rho=r/r_max in [0, 1]

            Returns:
                    (float): Wavefront error
    '''

    summation = 0
    for el in indices_coeffs:
        j = el[0]
        coeff = el[1]
        summation += coeff * z(j, rho, theta, rho_max)

    return summation


def get_coeff_from_rms(rms, indices):
    '''
    Returns the coefficient z_j = z of a sum of Zernike polynomials assuming the same coefficient for each contribution.

            Parameters:
                    rms (float): Total root mean square wavefront aberration
                    indices (list): List of tuples of the kind (j,), where j is the Zernike poylnomial index according to Noll's convention

            Returns:
                    (float): Zernike polynomial coefficient
    '''

    return abs(rms/np.sqrt(len(indices)))


def get_distribution_from_rms(rms, indices, slope):
    '''
    Returns the coefficient z_n^m = z of a sum of Zernike polynomials assuming the same coefficient for each contribution.

            Parameters:
                    rms (float): Total root mean square wavefront aberration
                    indices (list): List of tuples of the kind (j,), where j is the Zernike poylnomial index according to Noll's convention, should be ordered from lowest to highest
                    slope (float): slope of the distribution of Zernike polynomial coefficients, i. e. difference between coefficient z_i and z_(i+1)

            Returns:
                    (list): List of tuples of the kind (j, coefficient)
    '''
    
    # define function to be solved numerically to find first coefficient z_1
    def func(z_1):
        summation = 0
        for i in range(len(indices)):
            summation += (z_1 - i * slope)**2
            
        return np.sqrt(summation) - rms

    # solve for z_1
    z_1 = fsolve(func, 1)[0]
    
    # define all coefficients by iteratively subtracting the slope from the previouse coefficient, starting with z_1
    for counter, el in enumerate(indices):
        coeff = z_1 - counter * slope
            
        indices[counter] += (coeff,)

    return indices