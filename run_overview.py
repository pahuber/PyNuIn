#%%
import numpy as np
#import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe, get_coeff_from_rms, r, noll_index, get_distribution_from_rms
from pynuin.main.optics import aperture
from pynuin.util.plot import plot_zernike_wfe, plot_zernike_distribution
from numpy.fft import fft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker

# list_wfe = [(1, -1, 1), (1, 1, 1), (2, 0, 1), (2, 2, 1), (2, -2, 1), (3, 0, 1), (3, 3, 1), (3, -3, 1), (3, 1, 1), (3, -1, 1), (4, 0, 1), (4, 2, 1), (4, -2, 1), (4, 4, 1), (4, -4, 1)]

list_wfe = [(4,), (5,), (6,)]
rms = 1





list_wfe_new = get_distribution_from_rms(rms, list_wfe, 3)

plot_zernike_wfe(list_wfe_new)
plot_zernike_distribution(list_wfe_new)






# coeff = get_coeff_from_rms(1, list_wfe)
# for index in range(len(list_wfe)):
#     list_wfe[index] += (coeff,)

# plot_zernike_wfe(list_wfe)
# plot_zernike_distribution(list_wfe)





