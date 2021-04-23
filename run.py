import numpy as np
import matplotlib.pyplot as plt
from zernike import RZern
import poppy as pp
from matplotlib import cm
from mpl_toolkits import mplot3d


m = 1
n = 1

npix = 1000


def z(n, m, npix=1000):
    return pp.zernike.zernike(n, m, npix)


z11 = z(1, 1)
z22 = z(2, 2)


# 2D

# plt.imshow(z(1, 1) + z(2, 2) + z(2, 0))
# plt.colorbar()
# plt.show()


# 3D
# fig = plt.figure()
# ax = plt.axes(projection='3d')

# ax.plot_surface(length, length, z,cmap='viridis', edgecolor='none')
# ax.set_title('Surface plot')
# plt.show()



radius = 1
fov = 1
wavelength = 460e-9
pixscale = 0.01


osys = pp.OpticalSystem()
circular_aperture = pp.CircularAperture(radius=radius)
osys.add_pupil(circular_aperture)
zwfe = pp.ZernikeWFE(coefficients=[0, 0, 0, 0, 0], radius=radius)
print(zwfe)
osys.add_pupil(zwfe)
osys.add_detector(pixelscale=pixscale, fov_arcsec=fov)

psf = osys.calc_psf(wavelength=wavelength, display=True)[0].data

print(psf[100][100])


# plt.imshow(psf)
# plt.colorbar()
# plt.show()