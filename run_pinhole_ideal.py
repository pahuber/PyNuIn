#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe
from pynuin.main.optics import aperture
from numpy.fft import fft2, ifft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker


'''definitions'''
lam = 1 #1e-5 #m
D1 = 2 #0.01 #m
D2 = 2 * 1.22 * (lam/D1) #lam/D (units)

D1_2 = D1/2 #m
D2_2 = D2/2 #lam/D (units)


# aperturex = aperture(D1)


'''prepare matrices'''
xymin = -5
xymax = 5
x = np.linspace(xymin, xymax, 100)#np.arange(-7, 7, 0.2)
y = np.linspace(xymin, xymax, 100)#np.arange(-7, 7, 0.2)

wfe_img = np.zeros((len(x), len(y)))
aperture1 = np.zeros((len(x), len(y)))
aperture1norm = np.zeros((len(x), len(y)))
aperture2prime = np.zeros((len(x), len(y)))
e_ideal = np.zeros((len(x), len(y)), dtype=complex)
e_aberrated = np.zeros((len(x), len(y)), dtype=complex)
e_sum = np.zeros((len(x), len(y)), dtype=complex)
e_diff = np.zeros((len(x), len(y)), dtype=complex)


'''amplitude normalization'''
for counterx, elx in enumerate(x):
    for countery, ely in enumerate(y):
        
        # perform transformation to polar coordinates
        ra = np.sqrt(elx**2+ely**2)
        
        # define aprture 1
        if ra <= D1_2:
            aperture1norm[counterx][countery] = 1
            
            
            
a0 = abs(np.sqrt(1/np.sum(aperture1norm))) # such that the intensity over the aperture 1 is unity
# a0 = 1/(np.pi*D1_2**2) falsch


# print(np.sum(aperture1norm/len(x)))
# print(np.pi*D1_2**2)

'''calculations'''
for counterx, elx in enumerate(x):
    for countery, ely in enumerate(y):
        
        # perform transformation to polar coordinates
        ra = np.sqrt(elx**2+ely**2)
        the = np.arctan2(ely, elx)
        
        # define aprture 1
        if ra <= D1_2:
            aperture1[counterx][countery] = a0
            
        # define aperture 2
        if ra <= D2_2:
            aperture2prime[counterx][countery] = 1
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
        
# print initial intensity
print(np.sum(aperture1**2))

'''plotting'''
fig, axs = plt.subplots(3, 2)
extent = [xymin, xymax, xymin, xymax]

# aperture 1
img1 = axs[0, 0].imshow(aperture1, extent=extent)
fig.colorbar(img1, ax=axs[0, 0], fraction=0.046, pad=0.04)
axs[0, 0].set_title("A$_1$")

# aperture 2
img2 = axs[0, 1].imshow(aperture2prime, extent=extent)
fig.colorbar(img2, ax=axs[0, 1], fraction=0.046, pad=0.04)
axs[0, 1].set_title("A$_2$")

# e field at aperture 1
img3 = axs[1, 0].imshow(e_field_a1.real, extent=extent)
# img3.set_clim(1e-5, np.amax(intensity_a1))
fig.colorbar(img3, ax=axs[1, 0], fraction=0.046, pad=0.04)
axs[1, 0].set_title("$E_1$ Real")#"$\mathcal{R}(\mathcal{F}\{A_1\})$")

# e field at aperture 2
img4 = axs[1, 1].imshow(e_field_a2.real, extent=extent)
# img4.set_clim(1e-5, np.amax(intensity_a2))
fig.colorbar(img4, ax=axs[1, 1], fraction=0.046, pad=0.04)
axs[1, 1].set_title("$E_2$ Real") #$\mathcal{R}(\mathcal{F}\{A_1\}\cdot A_2)$")

# difference of e fields, i. e. what got filtered
img6 = axs[2, 0].imshow((e_field_a1 - e_field_a2).real, extent=extent)
# img4.set_clim(0.5e1, np.amax(intensity))
fig.colorbar(img6, ax=axs[2, 0], fraction=0.046, pad=0.04)
axs[2, 0].set_title("$E_1 - E_2$")#"$\mathcal{R}(\Delta E)$")

# intensity in image plane
img5 = axs[2, 1].imshow(intensity, extent=extent)
# img4.set_clim(0.5e1, np.amax(intensity))
fig.colorbar(img5, ax=axs[2, 1], fraction=0.046, pad=0.04)
intensity_final = round(np.sum(intensity), 2)
axs[2, 1].set_title("$I_{fin} = $" + str(intensity_final) + " $I_{init}$")#"$|\mathcal{R}(\mathcal{F}^{-1}\{\mathcal{F}\{A_1\}\cdot A_2\})|^2$")

# print final intensity as a fraction of the initial intensity
print(np.sum(intensity))

# plot
# plt.tight_layout(pad=1.5)
plt.subplots_adjust(wspace=-0.2, hspace=0.71)
plt.savefig("plot.pdf")
plt.show()