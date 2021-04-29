import numpy as np
import matplotlib.pyplot as plt
from zernike import RZern
import poppy as pp
from matplotlib import cm
from mpl_toolkits import mplot3d
import sympy as sym
from math import factorial
import pylab as pb

from main.model import R, Z, WFE, A_id, A_ab



# prepare matrices
x = np.arange(-1, 1, 0.05)
y = np.arange(-1, 1, 0.05)
amp_plus = np.zeros((len(x), len(y)))
amp_minus = np.zeros((len(x), len(y)))
i_max = np.zeros((len(x), len(y)))
i_min = np.zeros((len(x), len(y)))
null = np.zeros((len(x), len(y)))





# calculate matrices
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