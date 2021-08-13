#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe, get_coeff_from_rms, r, noll_index, get_distribution_from_rms
from pynuin.main.optics import aperture
from pynuin.util.plot import plot_zernike_wfe, plot_zernike_distribution, plot_aberrated_psf
from numpy.fft import fft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker


for num in range(1, 22):
    list_wfe = [(num, 0.1)]
    
    plot_zernike_wfe(list_wfe, name_plain="output/wfe_psf/wfe_" + str(num) + ".pdf")
    
    plot_aberrated_psf(list_wfe, name_plain="output/wfe_psf/psf_" + str(num) + ".pdf")




