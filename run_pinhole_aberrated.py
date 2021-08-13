#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe, get_coeff_from_rms
from pynuin.main.optics import aperture
from numpy.fft import fft2, ifft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker
from scipy.optimize import fsolve


'''definitions'''
lam = 1 #1e-5 #m
D1 = 0.5 #0.01 #m
D2 = 0.2 * (lam/D1) #lam/D (units)

# list_wfe = [(2, -2), (2, 2)]

# rms = lam/10


list_wfe = []
for i in range(5, 21+1, 1):
    list_wfe.append((i,))

rms = 0.1



D1_2 = D1/2 #m
D2_2 = D2/2 #lam/D (units)





'''rms to zernike coefficient calculations'''
coeff = get_coeff_from_rms(rms, list_wfe)
for index in range(len(list_wfe)):
    list_wfe[index] += (coeff,)
    

# aperturex = aperture(D1, list_wfe=list_wfe)

'''prepare matrices'''
xymin = -5
xymax = 5
x = np.linspace(xymin, xymax, 100)#np.arange(-7, 7, 0.2)
y = np.linspace(xymin, xymax, 100)#np.arange(-7, 7, 0.2)

wfe_img = np.zeros((len(x), len(y)))
aperture1_plot = np.zeros((len(x), len(y)))
aperture1 = np.zeros((len(x), len(y)), dtype=complex)
aperture1norm = np.zeros((len(x), len(y)), dtype=complex)
aperture2prime = np.zeros((len(x), len(y)))
e_ideal = np.zeros((len(x), len(y)), dtype=complex)
e_aberrated = np.zeros((len(x), len(y)), dtype=complex)
e_sum = np.zeros((len(x), len(y)), dtype=complex)
e_diff = np.zeros((len(x), len(y)), dtype=complex)



# '''amplitude normalization'''
# for counterx, elx in enumerate(x):
#     for countery, ely in enumerate(y):
        
#         # perform transformation to polar coordinates
#         ra = np.sqrt(elx**2+ely**2)
#         the = np.arctan2(ely, elx)
        
#         # specify wavefront error
#         wfe_gen = float(wfe(list_wfe, ra, the, D1_2))
        
#         # define aprture 1
#         aperture1norm[counterx][countery] = 1*np.heaviside(D1_2-ra, 1)*np.heaviside(ra, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)

# # normalize this amplitude to unit intensity
# amplitude_temp = np.sum(aperture1norm)

# # choose initial amplitude imaginary part because of degeneracy
# a0_imag = amplitude_temp.imag

# # define function to find roots for
# def func(a0_real):
    
#     func = 0
    
#     for el1 in aperture1norm:
#         for el2 in el1:
#             z_real = el2.real
#             z_imag = el2.imag
            
#             func += abs(a0_real*z_real - a0_imag*z_imag)**2
            
#     return func - 1

# # solve for roots, i. e. real part of amplitude a0
# a0_real = fsolve(func, 1)[0]

# # define full complex amplitude a0
# a0 = a0_real + a0_imag*1j


a0 = 1


'''calculations'''
for counterx, elx in enumerate(x):
    for countery, ely in enumerate(y):
        
        # perform transformation to polar coordinates
        ra = np.sqrt(elx**2+ely**2)
        the = np.arctan2(ely, elx)
        
        # # specify wafefront error
        wfe_gen = float(wfe(list_wfe, ra, the, D1_2))
        wfe_img[counterx][countery] = float(wfe(list_wfe, ra, the, D1_2))
        
        # define aprture 1
        aperture1_plot[counterx][countery] = 1*np.heaviside(D1_2-ra, 1)*np.heaviside(ra, 1)
        aperture1[counterx][countery] = a0*np.heaviside(D1_2-ra, 1)*np.heaviside(ra, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)
        # if ra <= D1_2:
        #     aperture1[counterx][countery] = a0
        
        # define aperture 2
        if ra <= D2_2:
            aperture2prime[counterx][countery] = 1
        
        # e field in aperture 1 plane
        e_field_a1 = fftshift(fft2(aperture1))
        intensity_a1 = abs(e_field_a1)**2
        
        # e field in aperture 2 plane
        e_field_a2 = e_field_a1 * aperture2prime
        intensity_a2 = abs(e_field_a2)**2
        
        # e field in imaging plane
        e_field = ifft2(e_field_a2)
        
        # calculate intensity in image plane
        intensity = abs(e_field)**2
        

        # # ideal e field
        # e_ideal[counterx][countery] = a0*np.heaviside(D_2-ra, 1)*np.heaviside(ra, 1)
        
        # # aberrated e field
        # e_aberrated[counterx][countery] = a0*np.heaviside(D_2-ra, 1)*np.heaviside(ra, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)

# define e fields and irradiances
# e_sum = e_ideal + e_aberrated
# e_diff = e_ideal - e_aberrated
# irr_ideal = abs(fftshift(fft2(e_ideal)).real)**2
# irr_aberrated = abs(fftshift(fft2(e_aberrated)).real)**2
# irr_max = abs(fftshift(fft2(e_sum)).real)**2
# irr_min = abs(fftshift(fft2(e_diff)).real)**2

# define null
# null = irr_min/irr_max
        
# print initial intensity
print(np.sum(abs(aperture1)**2))


'''plotting'''
fig, axs = plt.subplots(3, 2)
extent = [xymin, xymax, xymin, xymax]

# ideal irradiance
img1 = axs[0, 0].imshow(aperture1.real, extent=extent)
fig.colorbar(img1, ax=axs[0, 0], fraction=0.046, pad=0.04)
axs[0, 0].set_title("A$_1$")

# maximum irradiance
img2 = axs[0, 1].imshow(aperture2prime, extent=extent)
fig.colorbar(img2, ax=axs[0, 1], fraction=0.046, pad=0.04)
axs[0, 1].set_title("A$_2$")

# wfe
img3 = axs[1, 0].imshow(e_field_a1.real, extent=extent)
# img3.set_clim(1e-5, np.amax(intensity_a1))
fig.colorbar(img3, ax=axs[1, 0], fraction=0.046, pad=0.04)
axs[1, 0].set_title("$E_1$ Real")#"$\mathcal{R}(\mathcal{F}\{A_1\})$")


# wfe
img4 = axs[1, 1].imshow(e_field_a2.real, extent=extent)
# img4.set_clim(1e-5, np.amax(intensity_a2))
fig.colorbar(img4, ax=axs[1, 1], fraction=0.046, pad=0.04)
axs[1, 1].set_title("$E_2$ Real") #$\mathcal{R}(\mathcal{F}\{A_1\}\cdot A_2)$")



# difference between intensities
img6 = axs[2, 0].imshow((e_field_a1 - e_field_a2).real, extent=extent)
# img4.set_clim(0.5e1, np.amax(intensity))
fig.colorbar(img6, ax=axs[2, 0], fraction=0.046, pad=0.04)
axs[2, 0].set_title("$E_1 - E_2$")#"$\mathcal{R}(\Delta E)$")


# difference between intensities
img5 = axs[2, 1].imshow(intensity, extent=extent)
# img4.set_clim(0.5e1, np.amax(intensity))
fig.colorbar(img5, ax=axs[2, 1], fraction=0.046, pad=0.04)
intensity_final = round(np.sum(intensity), 2)
axs[2, 1].set_title("$I_{fin} = $" + str(intensity_final) + " $I_{init}$")#"$|\mathcal{R}(\mathcal{F}^{-1}\{\mathcal{F}\{A_1\}\cdot A_2\})|^2$")

# print final intensity as a fraction of the initial intensity
print(np.sum(intensity))

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

# plt.tight_layout(pad=1.5)
plt.subplots_adjust(wspace=-0.2, hspace=0.71)
plt.savefig("plot.pdf")
plt.show()


toplot = intensity_a1

plt.imshow(toplot)
plt.imsave("output/1.pdf", toplot)


# plt.imshow(null, norm=LogNorm())
# plt.clim(-1e-18, 0)
# plt.colorbar()

# plt.show()
# # cb = fig.colorbar(img7, ax=axs[1, 2], fraction=0.046, pad=0.04)
# # tick_locator = ticker.MaxNLocator(nbins=2)
# # cb.locator = tick_locator
# # cb.update_ticks()
# # axs[1, 2].set_title("Null")