import numpy as np
import sympy as sym
from math import factorial
from scipy.integrate import dblquad
from scipy.special import j1

        
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


# method for simple circular aperture function, returns 1 if rho in [0, 1), 0 else
# def circular_aperture(rho, theta=sym.Symbol("theta")):
#     return np.heaviside(1-rho, 1)*np.heaviside(rho, 1)




def integrand_id(rho,
         theta,
         e0, 
         t, 
         x, 
         y, 
         z,
         D, 
         lam, 
         c):
    
    
    R = np.sqrt(x**2+y**2+z**2)
    
    
    # print("rho " + str(rho))
    # print("theta " + str(theta))
    # print("e0 " + str(e0))
    # print("t " + str(t))
    # print("x " + str(x))
    # print("y " + str(y))
    # print("z " + str(z))
    # print("D " + str(D))
    # print("lam " + str(lam))
    # print("\n")
    
    
    aperture_fct = np.heaviside(D/2-rho, 1)*np.heaviside(rho, 1)
    
    term = np.exp(2*np.pi*1j*rho/(lam*R)*(x*np.cos(theta)+y*np.sin(theta)))
    
    # print(term)
    # print(aperture_fct)
    
    return aperture_fct*term*rho


# method to calculate E-field after passing through a circular aperture
def e_id(e0,
         t, 
         x, 
         y, 
         z,  
         D=1e-3, 
         lam=450e-9,
         rho_max=10*1e-3, #10*D
         c=3e8):
    
    R = np.sqrt(x**2+y**2+z**2)
    
    pre_term = e0/R*np.exp(2*np.pi*1j/lam*(c*t-R))
    
    integral = dblquad(integrand_id, 0, np.pi, lambda rho: 0, lambda rho: rho_max, args=(e0, t, x, y, z, D, lam, c))
    
    integral2 = dblquad(integrand_id, np.pi, 2*np.pi, lambda rho: 0, lambda rho: rho_max, args=(e0, t, x, y, z, D, lam, c))
    
    # integral = 0
    
    # steps = 2000
    
    # for rho in np.linspace(0, rho_max, steps):
        
    #     rho += rho_max/(steps*2)
    #     # print(rho)
    #     for theta in np.linspace(0, 2*np.pi, steps):
    #         theta += 2*np.pi/(steps*2)
    #         # print(rho)
    #         # print(theta)
    #         # print("\n")
            
    #         integral += integrand_id(rho, theta, e0, t, x, y, z, D, lam, c)

    
    e_field = pre_term * (integral[0] + integral2[0])
    print(integral[0])
    print(integral[1])
    
    print(integral2[0])
    print(integral2[1])
    
    return e_field


# method for analytical calculation of ideal e-field with bessel functions
def e_id_polar(t, 
         q, 
         z,
         e0,
         D=1e-2, 
         lam=450e-9,
         c=3e8):
    
    R = np.sqrt(q**2+z**2)
    
    pre_term = e0*lam*D/(2*q)*np.exp(2*np.pi*1j/lam*(c*t-R))
    
    bessel = j1(np.pi*q*D/(lam*R))

    # pre_term = 1
    
    e_field = pre_term * bessel
    
    return e_field






