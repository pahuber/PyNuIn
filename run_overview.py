#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe, get_coeff_from_rms, r, noll_index, get_distribution_from_rms
from pynuin.main.optics import aperture
from pynuin.util.plot import plot_zernike_wfe, plot_zernike_distribution, plot_aberrated_psf
from numpy.fft import fft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker

# list_wfe = [(1, -1, 1), (1, 1, 1), (2, 0, 1), (2, 2, 1), (2, -2, 1), (3, 0, 1), (3, 3, 1), (3, -3, 1), (3, 1, 1), (3, -1, 1), (4, 0, 1), (4, 2, 1), (4, -2, 1), (4, 4, 1), (4, -4, 1)]

# list_wfe = [(5,), (6,), (7,)]

list_wfe = []
for i in range(5, 21+1, 1):
    list_wfe.append((i,))

rms = 0.1


# list_wfe = [(4, 0.1)]

# unequal distribution with slope
# list_wfe = get_distribution_from_rms(rms, list_wfe, rms/23.600725)

# coeff_10 = list_wfe[0][1]
# coeff_21 = list_wfe[-1][1]

# print(coeff_10)
# print(coeff_21)

# print(coeff_10/coeff_21)

# print(list_wfe)

# check rms
# rms_check = 0
# for el in list_wfe:
#     coeff = el[1]
#     rms_check += coeff**2    
# rms_check = np.sqrt(rms_check)
# print(rms_check)







# equally distributed
coeff = get_coeff_from_rms(0.1, list_wfe)
for index in range(len(list_wfe)):
    list_wfe[index] += (coeff,)



plot_zernike_wfe(list_wfe, name_plain="output/wfe_psf/wfe_example.pdf")
    
plot_aberrated_psf(list_wfe, name_plain="output/wfe_psf/psf_example.pdf")


# plot_zernike_wfe(list_wfe)
# plot_zernike_distribution(list_wfe, normalized=True)
    
# plot_aberrated_psf(list_wfe, name_plain="z1.pdf")




