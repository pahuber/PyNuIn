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






# declare symbolic variables
t, lam, rho, theta = sym.symbols("t lam rho theta")



# wfe1 = WFE(piston=35e-9, tiltx=35e-9)


# a_id = A_id(10, t=10, lam=450e-9)
# a_ab = A_ab(10, wfe1, t=10, lam=450e-9)

# a_plus = a_id + a_ab

# a_minus = a_id - a_ab


# i_max = abs(a_plus)**2


# # print(a_plus)

# i_min = abs(a_minus)**2

# n = i_min/i_max




# prepare image
x = np.arange(-1, 1, 0.05)
y = np.arange(-1, 1, 0.05)
amp_plus = np.zeros((len(x), len(y)))
amp_minus = np.zeros((len(x), len(y)))
i_max = np.zeros((len(x), len(y)))
i_min = np.zeros((len(x), len(y)))
null = np.zeros((len(x), len(y)))

# calculate matrix
for counterx, elx in enumerate(x):
    for countery, ely in enumerate(y):
                
        # perform transformation to polar coordinates
        ra = np.sqrt(elx**2+ely**2)
        the = abs(np.arctan(ely/elx))
        
        wfe1 = float(WFE(piston=1e-8, tiltx=1e-8, tilty=1e-8, defocus=1e-9, rho=ra, theta=the))
        a_id = A_id(1, t=0, lam=450e-9)
        a_ab = A_ab(1, wfe1, t=0, lam=450e-9, rho=ra, theta=the)
        
        # print("-------------------")
        # print(a_ab)
        # print("-------------------")
        
        
        
        amp_plus[counterx][countery] = a_id + a_ab
        amp_minus[counterx][countery] = a_id - a_ab 
        i_max[counterx][countery] = abs(a_id + a_ab)**2
        i_min[counterx][countery] = abs(a_id - a_ab)**2
        null[counterx][countery] = (abs(a_id - a_ab)**2)/(abs(a_id + a_ab)**2)
        
        # img[counterx][countery] = z(3, -1, ra, the)
        
plt.imshow(null)
# plt.imshow(amp_minus)
       
plt.colorbar()
plt.show()

##################################





# ###################################
# #plot zernikes

# x = np.arange(-1, 1, 0.05)
# y = np.arange(-1, 1, 0.05)

# img = np.zeros((len(x), len(y)))


# for counterx, elx in enumerate(x):
#     for countery, ely in enumerate(y):
                
#         ra = np.sqrt(elx**2+ely**2)
#         the = abs(np.arctan(ely/elx))
        
#         print(countery)
        
#         img[counterx][countery] = WFE(rho=ra, theta=the, piston=1, tiltx=1, tilty=0, astig_oblique=0, defocus=0, astig_vertical=0)
#         # img[counterx][countery] = z(3, -1, ra, the)
        
# plt.imshow(img)
       
# plt.show()

# ##################################







# m = 1
# n = 1

# npix = 1000


# def z(n, m, npix=1000):
#     return pp.zernike.zernike(n, m, npix)


# z11 = z(1, 1)
# z22 = z(2, 2)


# # 2D

# plt.imshow(z(1, 1) + z(2, 2) + z(2, 0))
# plt.colorbar()
# plt.show()


# # 3D
# # fig = plt.figure()
# # ax = plt.axes(projection='3d')

# # ax.plot_surface(length, length, z,cmap='viridis', edgecolor='none')
# # ax.set_title('Surface plot')
# # plt.show()



# radius = 1
# fov = 1
# wavelength = 460e-9
# pixscale = 0.01


# osys = pp.OpticalSystem()
# circular_aperture = pp.CircularAperture(radius=radius)
# osys.add_pupil(circular_aperture)
# zwfe = pp.ZernikeWFE(coefficients=[0, 0, 0, 0, 0], radius=radius)
# print(zwfe)
# osys.add_pupil(zwfe)
# osys.add_detector(pixelscale=pixscale, fov_arcsec=fov)

# psf = osys.calc_psf(wavelength=wavelength, display=True)[0].data

# # print(psf[100][100])


# # plt.imshow(psf)
# # plt.colorbar()
# # plt.show()