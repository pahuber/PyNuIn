#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.model import wfe
from numpy.fft import fft2, ifft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker


'''prepare matrices'''
x = np.linspace(-5, 5, 100)#np.arange(-7, 7, 0.2)
y = np.linspace(-5, 5, 100)#np.arange(-7, 7, 0.2)

wfe_img = np.zeros((len(x), len(y)))
aperture1 = np.zeros((len(x), len(y)))
aperture2prime = np.zeros((len(x), len(y)))
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
        D1_2 = 0.8
        
        D2_2 = 0.8
        lam = 1e-5
        # rho_max = np.sqrt(max(x)**2+max(y)**2)
        
        
        # define aprture 1
        if ra <= D1_2:
            aperture1[counterx][countery] = a0
            
        # define aperture 2
        if ra <= D2_2:
            aperture2prime[counterx][countery] = a0
        aperture2 = fftshift(fft2(aperture2prime))
        # aperture2 = abs(aperture2.real)**2
        
        # e field in aperture 1 plane
        e_field_a1 = fftshift(fft2(aperture1))
        intensity_a1 = abs(e_field_a1.real)**2
        
        # e field in aperture 2 plane
        e_field_a2 = e_field_a1 * aperture2prime
        intensity_a2 = abs(e_field_a2.real)**2
        
        # e field in imaging plane
        e_field = ifft2(e_field_a2)
        
        # calculate intensity in image plane
        intensity = abs(e_field.real)**2
        

'''plotting'''
fig, axs = plt.subplots(3, 2)

# aperture 1
img1 = axs[0, 0].imshow(aperture1)
fig.colorbar(img1, ax=axs[0, 0], fraction=0.046, pad=0.04)
axs[0, 0].set_title("A$_1$")

# aperture 2
img2 = axs[0, 1].imshow(aperture2prime)
fig.colorbar(img2, ax=axs[0, 1], fraction=0.046, pad=0.04)
axs[0, 1].set_title("A$_2$")

# e field at aperture 1
img3 = axs[1, 0].imshow(e_field_a1.real)
# img3.set_clim(1e-5, np.amax(intensity_a1))
fig.colorbar(img3, ax=axs[1, 0], fraction=0.046, pad=0.04)
axs[1, 0].set_title("$\mathcal{R}(\mathcal{F}\{A_1\})$")

# e field at aperture 2
img4 = axs[1, 1].imshow(e_field_a2.real)
# img4.set_clim(1e-5, np.amax(intensity_a2))
fig.colorbar(img4, ax=axs[1, 1], fraction=0.046, pad=0.04)
axs[1, 1].set_title("$\mathcal{R}(\mathcal{F}\{A_1\}\cdot A_2)$")

# difference of e fields, i. e. what got filtered
img6 = axs[2, 0].imshow((e_field_a1 - e_field_a2).real)
# img4.set_clim(0.5e1, np.amax(intensity))
fig.colorbar(img6, ax=axs[2, 0], fraction=0.046, pad=0.04)
axs[2, 0].set_title("$\mathcal{R}(\Delta E)$")

# intensity in image plane
img5 = axs[2, 1].imshow(intensity)
# img4.set_clim(0.5e1, np.amax(intensity))
fig.colorbar(img5, ax=axs[2, 1], fraction=0.046, pad=0.04)
axs[2, 1].set_title("$|\mathcal{R}(\mathcal{F}^{-1}\{\mathcal{F}\{A_1\}\cdot A_2\})|^2$")

# plot
# plt.tight_layout(pad=1.5)
plt.subplots_adjust(wspace=-0.2, hspace=0.71)
plt.savefig("plot.pdf")
plt.show()