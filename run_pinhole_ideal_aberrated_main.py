#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe, get_coeff_from_rms
from pynuin.main.optics import aperture, normalize_apertures, normalize_apertures2
from numpy.fft import fft2, ifft2, fftshift
# from scipy.fftpack import fft2, fftshift, ifft2
from matplotlib.colors import LogNorm
from matplotlib import ticker
from scipy.optimize import fsolve


'''computational specifications'''
n = 100
size = n**2


'''physical specificaations'''
lam = 1 #1e-5 #m
D1 = 20 #0.01 #m
D2 = 2 * 1.22 * (lam/D1) * n #lam/D (units)
list_wfe = [(2, 0), (2, 2), (2, -2), (3, 1), (3, -1), (3, 3), (3, -3), (4, 0), (4, 2), (4, -2), (4, 4), (4, -4), (5, 5), (5, 3), (5, 1), (5, -5), (5, -3), (5, -1)]
rms = 0.021*lam

D1_2 = D1/2 #m
D2_2 = D2/2 #lam/D (units)


'''rms to zernike coefficient calculations'''
coeff = get_coeff_from_rms(rms, list_wfe)
for index in range(len(list_wfe)):
    list_wfe[index] += (coeff,)


'''create apertures'''
# create ideal aperture
a1_id = aperture(D = D1, 
                 lam = lam,
                 a0 = 1,
                 list_wfe = None,
                 n = n)

# print(np.sum(a1_id)/np.count_nonzero(a1_id))

# create aberrated aperture
a1_ab = aperture(D = D1, 
                 lam = lam,
                 a0 = 1,
                 list_wfe = list_wfe,
                 n = n)

# print(np.count_nonzero(a1_ab))
# print(np.sum(a1_id) + np.sum(a1_ab))
# print(np.sum(a1_id)/size/100)
# a1_id = a1_id/size/100
# a1_ab = a1_ab/size/100

# normalize apertures
# a1_id, a1_ab, a0 = normalize_apertures(a1_id, a1_ab, 1)
# a1_id, a1_ab, a0 = normalize_apertures2(a1_id, a1_ab, 1)

# print(a0)

# print(type(a0))

# print(((np.sum(a1_id) + np.sum(a1_ab)).real)**2)
# a1_id = a1_id/size
# a1_ab = a1_ab/size
# print(a1_id.shape)
# print(np.sum(a1_id))
# print(np.sum(a1_ab))
# print(a0)

# create aperture in pinhole plane
a2 = aperture(D = D2, 
              lam = lam,
              a0 = 1,
              list_wfe = None,
              n = n)


'''calculations'''
# e field in pinhole plane
e_pinhole_id = fftshift(fft2(a1_id))
e_pinhole_ab = fftshift(fft2(a1_ab))

# e field in final plane
e_pinhole_filtered_id = e_pinhole_id * a2
e_pinhole_filtered_ab = e_pinhole_ab  * a2
e_final_id = ifft2(e_pinhole_filtered_id)
e_final_ab = ifft2(e_pinhole_filtered_ab)

# sum and diff of e fields
e_plus = e_final_id + e_final_ab
e_minus = e_final_id - e_final_ab

# calculate intensity in final plane
intensity_max = abs(e_plus.real)**2
intensity_min = abs(e_minus.real)**2

# define null
imax = np.sum(intensity_max)#/intensity_max.size
imin = np.sum(intensity_min)#/intensity_max.size
null = imin/imax

# calculate initial common intensity, should equal intensity_init from above, i. e. I_init total (|E_1 + E_2|^2)
iinit_common = np.sum(abs((a1_id + a1_ab).real)**2)#/intensity_max.size

# prints
print(imax)
print(iinit_common)

print("Common initial intensity: " + str(iinit_common))
print("Null: " + str(null))
print("Throughput: " + str(round(imax/iinit_common * 100, 1)) + " %")


'''plotting'''
fig, axs = plt.subplots(3, 2)
xymin = -n/2
xymax = n/2
extent = [xymin, xymax, xymin, xymax]

# ideal aperture 1
img1 = axs[0, 0].imshow(a1_id.real, extent=extent)
fig.colorbar(img1, ax=axs[0, 0], fraction=0.046, pad=0.04)
axs[0, 0].set_title("A$_1$")

# aberrated aperture 1
img2 = axs[0, 1].imshow(a2.real, extent=extent)
fig.colorbar(img2, ax=axs[0, 1], fraction=0.046, pad=0.04)
axs[0, 1].set_title("A$_2$")

# filtered ideal e field 
img3 = axs[1, 0].imshow(e_pinhole_filtered_id.real, extent=extent)
# img3.set_clim(1e-5, np.amax(intensity_a1))
fig.colorbar(img3, ax=axs[1, 0], fraction=0.046, pad=0.04)
axs[1, 0].set_title("E id clean")

# filtered aberrated e field
img4 = axs[1, 1].imshow(e_pinhole_filtered_ab.real, extent=extent)
# img4.set_clim(1e-5, np.amax(intensity_a2))
fig.colorbar(img4, ax=axs[1, 1], fraction=0.046, pad=0.04)
axs[1, 1].set_title("E ab clean")

# intensity max
img6 = axs[2, 0].imshow(intensity_max, extent=extent)
# img4.set_clim(0.5e1, np.amax(intensity))
fig.colorbar(img6, ax=axs[2, 0], fraction=0.046, pad=0.04)
axs[2, 0].set_title("I max")

# intensity min
img5 = axs[2, 1].imshow(intensity_min, extent=extent)
# img4.set_clim(0.5e1, np.amax(intensity))
fig.colorbar(img5, ax=axs[2, 1], fraction=0.046, pad=0.04)
axs[2, 1].set_title("I min")

# plt.tight_layout(pad=1.5)
plt.subplots_adjust(wspace=-0.2, hspace=0.71)
plt.savefig("plot.pdf")
plt.show()