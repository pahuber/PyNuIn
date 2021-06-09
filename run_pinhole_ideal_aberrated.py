#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe, get_coeff_from_rms
from numpy.fft import fft2, ifft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker
from scipy.optimize import fsolve


'''definitions'''
lam = 1 #1e-5 #m
D1 = 1 #0.01 #m
D2 = 1.0 * (lam/D1) #lam/D (units)
list_wfe = [(2, 0), (2, 2), (2, -2), (3, 1), (3, -1), (3, 3), (3, -3), (4, 0), (4, 2), (4, -2), (4, 4), (4, -4), (5, 5), (5, 3), (5, 1), (5, -5), (5, -3), (5, -1)]
rms = 0.05*lam

D1_2 = D1/2 #m
D2_2 = D2/2 #lam/D (units)


'''rms to zernike coefficient calculations'''
coeff = get_coeff_from_rms(rms, list_wfe)
for index in range(len(list_wfe)):
    list_wfe[index] += (coeff,)


'''prepare matrices'''
xymin = -5
xymax = 5
x = np.linspace(xymin, xymax, 300)#np.arange(-7, 7, 0.2)
y = np.linspace(xymin, xymax, 300)#np.arange(-7, 7, 0.2)

wfe_img = np.zeros((len(x), len(y)))
aperture1_id = np.zeros((len(x), len(y)), dtype=complex)
aperture1_ab = np.zeros((len(x), len(y)), dtype=complex)
aperture1_id_norm = np.zeros((len(x), len(y)), dtype=complex)
aperture1_ab_norm = np.zeros((len(x), len(y)), dtype=complex)
aperture1_id_norm_check = np.zeros((len(x), len(y)), dtype=complex)
aperture1_ab_norm_check = np.zeros((len(x), len(y)), dtype=complex)
aperture2 = np.zeros((len(x), len(y)))
e_ideal = np.zeros((len(x), len(y)), dtype=complex)
e_aberrated = np.zeros((len(x), len(y)), dtype=complex)
e_sum = np.zeros((len(x), len(y)), dtype=complex)
e_diff = np.zeros((len(x), len(y)), dtype=complex)


'''amplitude normalization common'''
for counterx, elx in enumerate(x):
    for countery, ely in enumerate(y):
        
        # perform transformation to polar coordinates
        ra = np.sqrt(elx**2+ely**2)
        the = np.arctan2(ely, elx)
        
        # specify wavefront error
        wfe_gen = float(wfe(list_wfe, ra, the, D1_2))
        
        # define aprture 1
        aperture1_id_norm[counterx][countery] = 1*np.heaviside(D1_2-ra, 1)*np.heaviside(ra, 1)
        aperture1_ab_norm[counterx][countery] = 1*np.heaviside(D1_2-ra, 1)*np.heaviside(ra, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)
        
# normalize this amplitude to unit intensity
aperture1_common  = aperture1_id_norm + aperture1_ab_norm
amplitude_temp = np.sum(aperture1_common)

# choose initial amplitude imaginary part because of degeneracy
a0_imag = amplitude_temp.imag

# define function to find roots for
def func(a0_real):
    
    func = 0
    
    for el1 in aperture1_common:
        for el2 in el1:
            z_real = el2.real
            z_imag = el2.imag
            
            func += abs(a0_real*z_real - a0_imag*z_imag)**2
            
    return func - 1

# solve for roots, i. e. real part of amplitude a0
a0_real = fsolve(func, 1)[0]

# define full complex amplitude a0
a0_common = a0_real + a0_imag*1j

print(a0_common)

# check initial intensity
# for counterx, elx in enumerate(x):
#     for countery, ely in enumerate(y):
        
#         # perform transformation to polar coordinates
#         ra = np.sqrt(elx**2+ely**2)
#         the = np.arctan2(ely, elx)
        
#         # specify wavefront error
#         wfe_gen = float(wfe(list_wfe, ra, the, D1_2))
            
aperture1_id_norm_check = a0_common*aperture1_id_norm
aperture1_ab_norm_check = a0_common*aperture1_ab_norm
        
iinit_common_check = np.sum(abs((aperture1_id_norm_check + aperture1_ab_norm_check).real)**2)
print(iinit_common_check)
if round(iinit_common_check, 3) != 1.0:
    print("Re-normalizing")
    
    
    # choose initial amplitude imaginary part because of degeneracy
    # a0_imag = amplitude_temp.imag
    
    # define function to find roots for
    def func(a0_imag_new):
        
        func = 0
        
        for el1 in aperture1_common:
            for el2 in el1:
                z_real = el2.real
                z_imag = el2.imag
                
                func += abs(a0_real*z_real - a0_imag_new*z_imag)**2
                
        return func - 1
    
    # solve for roots, i. e. real part of amplitude a0
    a0_imag_new = fsolve(func, 1)[0]
    
    # define full complex amplitude a0
    a0_common = a0_real + a0_imag_new*1j
    
    print(a0_common)
        



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
        # print(a0_common)
        aperture1_id[counterx][countery] = a0_common*np.heaviside(D1_2-ra, 1)*np.heaviside(ra, 1)
        aperture1_ab[counterx][countery] = a0_common*np.heaviside(D1_2-ra, 1)*np.heaviside(ra, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)
        # if ra <= D1_2:
        #     aperture1[counterx][countery] = a0
        
        # define aperture 2
        if ra <= D2_2:
            aperture2[counterx][countery] = 1

# e field in aperture 1 plane
e_field_a1_id = fftshift(fft2(aperture1_id))
e_field_a1_ab = fftshift(fft2(aperture1_ab))
# e_field_a1_plus = fftshift(fft2(aperture1_id+aperture1_ab)) #direct
# e_field_a1_minus = fftshift(fft2(aperture1_id-aperture1_ab)) #direct
# intensity_a1 = abs(e_field_a1.real)**2

# e field in aperture 2 plane
e_field_a2_id = e_field_a1_id * aperture2
e_field_a2_ab = e_field_a1_ab * aperture2
# e_field_a2_plus = e_field_a1_plus * aperture2 #direct
# e_field_a2_minus = e_field_a1_minus * aperture2 #direct
# intensity_a2 = abs(e_field_a2.real)**2

# e field in imaging plane
e_field_id = ifft2(e_field_a2_id)
e_field_ab = ifft2(e_field_a2_ab)
# e_field_dirty_id = ifft2(e_field_a1_id)
# e_field_dirty_ab = ifft2(e_field_a1_ab)

# e_field_plus = ifft2(e_field_a2_plus) #direct
# e_field_minus = ifft2(e_field_a2_minus) #direct

# sum and diff of e fields
e_plus = e_field_id + e_field_ab
e_minus = e_field_id - e_field_ab

# e_plus_dirty = e_field_dirty_id + e_field_dirty_ab
# e_minus_dirty = e_field_dirty_id - e_field_dirty_ab

# calculate intensity in image plane
intensity_max = abs(e_plus.real)**2
intensity_min = abs(e_minus.real)**2


# intensity_max_dirty = abs(e_plus_dirty.real)**2
# intensity_min_dirty = abs(e_minus_dirty.real)**2

# intensity_plus_direct = abs(e_field_plus.real)**2
# intensity_minus_direct = abs(e_field_minus.real)**2

# define null
imax = np.sum(intensity_max)
imin = np.sum(intensity_min)
null = imin/imax

print(imax)

# imaxd = np.sum(intensity_max_dirty)
# imind = np.sum(intensity_min_dirty)
# nulld = imind/imaxd
# # print(imax)
# print(imin)
# print(imax)

# I_init total (|E_1 + E_2|^2)
iinit_common = np.sum(abs((aperture1_id + aperture1_ab).real)**2)
print(iinit_common)

# I_init separat
# iinit = np.sum(abs(aperture1_id.real)**2) + np.sum(abs(aperture1_ab.real)**2)
# print(iinit)

# print(np.sum(abs(aperture1_id.real)**2) + np.sum(abs(aperture1_ab.real)**2)

# print(np.sum(abs(aperture1_id.real)**2)) #iinit id
# print(np.sum(abs(aperture1_ab.real)**2)) # iinit ab

# print(imax)
# print(iinit) #iinit sep

print("Common initial intensity: " + str(iinit_common))
print("Null: " + str(round(null, 10)))
print("Throughput: " + str(round(imax/iinit_common * 100, 1)) + " %")

# print(np.sum(intensity_minus_direct)/np.sum(intensity_plus_direct))
# null = irr_min/irr_max
        

'''plotting'''
fig, axs = plt.subplots(3, 2)
extent = [xymin, xymax, xymin, xymax]

# ideal irradiance
img1 = axs[0, 0].imshow(aperture1_id.real, extent=extent)
fig.colorbar(img1, ax=axs[0, 0], fraction=0.046, pad=0.04)
axs[0, 0].set_title("A$_1$")

# maximum irradiance
img2 = axs[0, 1].imshow(aperture2, extent=extent)
fig.colorbar(img2, ax=axs[0, 1], fraction=0.046, pad=0.04)
axs[0, 1].set_title("A$_2$")

# wfe
img3 = axs[1, 0].imshow(e_field_a2_id.real, extent=extent)
# img3.set_clim(1e-5, np.amax(intensity_a1))
fig.colorbar(img3, ax=axs[1, 0], fraction=0.046, pad=0.04)
axs[1, 0].set_title("E id clean")


# wfe
img4 = axs[1, 1].imshow(e_field_a2_ab.real, extent=extent)
# img4.set_clim(1e-5, np.amax(intensity_a2))
fig.colorbar(img4, ax=axs[1, 1], fraction=0.046, pad=0.04)
axs[1, 1].set_title("E ab clean")



# difference between intensities
img6 = axs[2, 0].imshow(intensity_max, extent=extent)
# img4.set_clim(0.5e1, np.amax(intensity))
fig.colorbar(img6, ax=axs[2, 0], fraction=0.046, pad=0.04)
axs[2, 0].set_title("I max")


# difference between intensities
img5 = axs[2, 1].imshow(intensity_min, extent=extent)
# img4.set_clim(0.5e1, np.amax(intensity))
fig.colorbar(img5, ax=axs[2, 1], fraction=0.046, pad=0.04)
axs[2, 1].set_title("I min")



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



# plt.imshow(null, norm=LogNorm())
# plt.clim(-1e-18, 0)
# plt.colorbar()

# plt.show()
# # cb = fig.colorbar(img7, ax=axs[1, 2], fraction=0.046, pad=0.04)
# # tick_locator = ticker.MaxNLocator(nbins=2)
# # cb.locator = tick_locator
# # cb.update_ticks()
# # axs[1, 2].set_title("Null")