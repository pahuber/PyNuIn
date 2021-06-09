#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe
from numpy.fft import fft2, ifft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker


'''prepare matrices'''
x = np.linspace(-5, 5, 100)#np.arange(-7, 7, 0.2)
y = np.linspace(-5, 5, 100)#np.arange(-7, 7, 0.2)

wfe_img = np.zeros((len(x), len(y)))
aperture1_plot = np.zeros((len(x), len(y)))
aperture1 = np.zeros((len(x), len(y)), dtype=complex)
aperture2prime = np.zeros((len(x), len(y)))
e_ideal = np.zeros((len(x), len(y)), dtype=complex)
e_aberrated = np.zeros((len(x), len(y)), dtype=complex)
e_sum = np.zeros((len(x), len(y)), dtype=complex)
e_diff = np.zeros((len(x), len(y)), dtype=complex)

imins = []
imaxs = []
nulls = []
apertures = np.arange(0.1, 2, 0.1)


'''calculations'''
for d2_2 in apertures:
    for counterx, elx in enumerate(x):
        for countery, ely in enumerate(y):
            
            # perform transformation to polar coordinates
            ra = np.sqrt(elx**2+ely**2)
            the = np.arctan2(ely, elx)
            
            # define constants
            a0 = 1
            D1_2 = 0.8
            
            D2_2 = d2_2
            lam = 1e-5
            
            # # specify wafefront error
            list_wfe = [(2, -2, 0.5e-5)]
            wfe_gen = float(wfe(list_wfe, ra, the))
            wfe_img[counterx][countery] = float(wfe(list_wfe, ra, the))
            
            # define aprture 1
            aperture1_plot[counterx][countery] = a0*np.heaviside(D1_2-ra, 1)*np.heaviside(ra, 1)
            aperture1[counterx][countery] = a0*np.heaviside(D1_2-ra, 1)*np.heaviside(ra, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)
            # if ra <= D1_2:
            #     aperture1[counterx][countery] = a0
            
            # define aperture 2
            if ra <= D2_2:
                aperture2prime[counterx][countery] = a0
            
            # e field in aperture 1 plane
            e_field_a1_id = fftshift(fft2(aperture1_plot))
            e_field_a1_ab = fftshift(fft2(aperture1))
            # intensity_a1 = abs(e_field_a1.real)**2
            
            # e field in aperture 2 plane
            e_field_a2_id = e_field_a1_id * aperture2prime
            e_field_a2_ab = e_field_a1_ab * aperture2prime
            # intensity_a2 = abs(e_field_a2.real)**2
            
            # e field in imaging plane
            e_field_id = ifft2(e_field_a2_id)
            e_field_ab = ifft2(e_field_a2_ab)
            
            # sum and diff of e fields
            e_plus = e_field_id + e_field_ab
            e_minus = e_field_id - e_field_ab
            
            # calculate intensity in image plane
            intensity_max = abs(e_plus.real)**2
            intensity_min = abs(e_minus.real)**2
            
    
    
    # define null
    imax = np.sum(intensity_max)
    imin = np.sum(intensity_min)
    null = imin/imax
    print(imax)
    print(imin)
    print(null)
    imins.append(imin)
    imaxs.append(imax)
    nulls.append(null)
# null = irr_min/irr_max
        


# plt.plot(apertures, imins, label="imins")
# plt.plot(apertures, imaxs, label="imaxs")
plt.plot(apertures, nulls, label="nulls")
plt.legend()
plt.show()


'''plotting'''
# fig, axs = plt.subplots(3, 2)

# # ideal irradiance
# img1 = axs[0, 0].imshow(aperture1_plot)
# fig.colorbar(img1, ax=axs[0, 0], fraction=0.046, pad=0.04)
# axs[0, 0].set_title("A$_1$")

# # maximum irradiance
# img2 = axs[0, 1].imshow(aperture2prime)
# fig.colorbar(img2, ax=axs[0, 1], fraction=0.046, pad=0.04)
# axs[0, 1].set_title("A$_2$")

# # wfe
# img3 = axs[1, 0].imshow(e_field_a2_id.real)
# # img3.set_clim(1e-5, np.amax(intensity_a1))
# fig.colorbar(img3, ax=axs[1, 0], fraction=0.046, pad=0.04)
# axs[1, 0].set_title("E id clean")


# # wfe
# img4 = axs[1, 1].imshow(e_field_a2_ab.real)
# # img4.set_clim(1e-5, np.amax(intensity_a2))
# fig.colorbar(img4, ax=axs[1, 1], fraction=0.046, pad=0.04)
# axs[1, 1].set_title("E ab clean")



# # difference between intensities
# img6 = axs[2, 0].imshow(intensity_max)
# # img4.set_clim(0.5e1, np.amax(intensity))
# fig.colorbar(img6, ax=axs[2, 0], fraction=0.046, pad=0.04)
# axs[2, 0].set_title("I max")


# # difference between intensities
# img5 = axs[2, 1].imshow(intensity_min)
# # img4.set_clim(0.5e1, np.amax(intensity))
# fig.colorbar(img5, ax=axs[2, 1], fraction=0.046, pad=0.04)
# axs[2, 1].set_title("I min")



# # irr aberrated
# img5 = axs[1, 0].imshow(irr_aberrated)
# fig.colorbar(img5, ax=axs[1, 0], fraction=0.046, pad=0.04)
# axs[1, 0].set_title("$I_{aberrated}$")

# # minimum irradiance
# img6 = axs[1, 1].imshow(irr_min)
# fig.colorbar(img6, ax=axs[1, 1], fraction=0.046, pad=0.04)
# axs[1, 1].set_title("$I_{min}$")

# # null
# # img7 = axs[1, 2].imshow(null, norm=LogNorm())
# # cb = fig.colorbar(img7, ax=axs[1, 2], fraction=0.046, pad=0.04)
# # tick_locator = ticker.MaxNLocator(nbins=2)
# # cb.locator = tick_locator
# # cb.update_ticks()
# axs[1, 2].set_title("Null")

# # plt.tight_layout(pad=1.5)
# plt.subplots_adjust(wspace=-0.2, hspace=0.71)
# plt.savefig("plot.pdf")
# plt.show()



# plt.imshow(null, norm=LogNorm())
# plt.clim(-1e-18, 0)
# plt.colorbar()

# plt.show()
# # cb = fig.colorbar(img7, ax=axs[1, 2], fraction=0.046, pad=0.04)
# # tick_locator = ticker.MaxNLocator(nbins=2)
# # cb.locator = tick_locator
# # cb.update_ticks()
# # axs[1, 2].set_title("Null")