#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.model import wfe
from numpy.fft import fft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker


'''prepare matrices'''
x = np.arange(-10, 10, 0.2)
y = np.arange(-10, 10, 0.2)

wfe_img = np.zeros((len(x), len(y)))
e_ideal = np.zeros((len(x), len(y)), dtype=complex)
e_aberrated = np.zeros((len(x), len(y)), dtype=complex)
e_sum = np.zeros((len(x), len(y)), dtype=complex)
e_diff = np.zeros((len(x), len(y)), dtype=complex)


'''calculations'''
for counterx, elx in enumerate(x):
    for countery, ely in enumerate(y):
        
        # perform transformation to polar coordinates
        ra = np.sqrt(elx**2+ely**2)
        the = np.arctan2(ely, elx)
        
        # define constants
        a0 = 1
        D_2 =0.5
        lam = 1e-5
        # rho_max = np.sqrt(max(x)**2+max(y)**2)
        
        # specify wafefront error
        list_wfe = [(1, 1, 0.5e-5), (1, -1, 0.5e-5)]
        wfe_gen = float(wfe(list_wfe, ra, the))
        wfe_img[counterx][countery] = float(wfe(list_wfe, ra, the))

        # ideal e field
        e_ideal[counterx][countery] = a0*np.heaviside(D_2-ra, 1)*np.heaviside(ra, 1)
        
        # aberrated e field
        e_aberrated[counterx][countery] = a0*np.heaviside(D_2-ra, 1)*np.heaviside(ra, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)

# define e fields and irradiances
e_sum = e_ideal + e_aberrated
e_diff = e_ideal - e_aberrated
irr_ideal = abs(fftshift(fft2(e_ideal)).real)**2
irr_aberrated = abs(fftshift(fft2(e_aberrated)).real)**2
irr_max = abs(fftshift(fft2(e_sum)).real)**2
irr_min = abs(fftshift(fft2(e_diff)).real)**2

# define null
null = irr_min/irr_max


'''plotting'''
fig, axs = plt.subplots(2, 3)

# ideal irradiance
img1 = axs[0, 0].imshow(irr_ideal)
fig.colorbar(img1, ax=axs[0, 0], fraction=0.046, pad=0.04)
axs[0, 0].set_title("$I_{ideal}$")

# maximum irradiance
img2 = axs[0, 1].imshow(irr_max)
fig.colorbar(img2, ax=axs[0, 1], fraction=0.046, pad=0.04)
axs[0, 1].set_title("$I_{max}$")

# wfe
img3 = axs[0, 2].imshow(wfe_img)
fig.colorbar(img3, ax=axs[0, 2], fraction=0.046, pad=0.04)
axs[0, 2].set_title("WFE")

# irr aberrated
img5 = axs[1, 0].imshow(irr_aberrated)
fig.colorbar(img5, ax=axs[1, 0], fraction=0.046, pad=0.04)
axs[1, 0].set_title("$I_{aberrated}$")

# minimum irradiance
img6 = axs[1, 1].imshow(irr_min)
fig.colorbar(img6, ax=axs[1, 1], fraction=0.046, pad=0.04)
axs[1, 1].set_title("$I_{min}$")

# null
img7 = axs[1, 2].imshow(null, norm=LogNorm())
cb = fig.colorbar(img7, ax=axs[1, 2], fraction=0.046, pad=0.04)
# tick_locator = ticker.MaxNLocator(nbins=2)
# cb.locator = tick_locator
# cb.update_ticks()
axs[1, 2].set_title("Null")

# plt.tight_layout(pad=1.5)
plt.subplots_adjust(wspace=0.6, hspace=0.21)
plt.savefig("plot.pdf")
plt.show()