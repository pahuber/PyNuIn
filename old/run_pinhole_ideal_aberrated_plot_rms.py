#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe, get_coeff_from_rms
from numpy.fft import fft2, ifft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker
from scipy.optimize import fsolve


'''definitions'''
D1_2 = 0.8
# D2_2 = 0.5
lam = 1e-5
# list_wfe = [(1, 1, 0.1e-5)]



pinholes = [0.2, 0.4]

rmss = [1e-7, 1e-6]

list_wfe = [(1, 1), (1, -1)]


# pinholes = [0.8, 1]

# rmss= [1e-5, 1e-4]

output_null = np.zeros((len(pinholes), len(rmss)))

output_throughput = np.zeros((len(pinholes), len(rmss)))





'''Loop through RMS's'''

for counter_rms, rms in enumerate(rmss):
    
    print("\n")
    print("----------------------------")
    print("RMS: " + str(rms))
    print("\n")
    
    
    '''prepare matrices for normalisation'''
    x = np.linspace(-5, 5, 100)#np.arange(-7, 7, 0.2)
    y = np.linspace(-5, 5, 100)#np.arange(-7, 7, 0.2)
    
    wfe_img = np.zeros((len(x), len(y)))
    aperture1_id = np.zeros((len(x), len(y)), dtype=complex)
    aperture1_ab = np.zeros((len(x), len(y)), dtype=complex)
    aperture1_id_norm = np.zeros((len(x), len(y)))
    aperture1_ab_norm = np.zeros((len(x), len(y)), dtype=complex)
    aperture2 = np.zeros((len(x), len(y)))
    e_ideal = np.zeros((len(x), len(y)), dtype=complex)
    e_aberrated = np.zeros((len(x), len(y)), dtype=complex)
    e_sum = np.zeros((len(x), len(y)), dtype=complex)
    e_diff = np.zeros((len(x), len(y)), dtype=complex)
    
    
    coeff = get_coeff_from_rms(rms, list_wfe)
    
    for index in range(len(list_wfe)):
        list_wfe[index] += (coeff,)
    
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
    
    
    
    '''Loop through pinholes'''
    
    for counter_pinhole, D2_2 in enumerate(pinholes):
    
        print("Pinhole Radius: " + str(D2_2))
        
        
    
        '''prepare matrices for calculations'''
        x = np.linspace(-5, 5, 100)#np.arange(-7, 7, 0.2)
        y = np.linspace(-5, 5, 100)#np.arange(-7, 7, 0.2)
        
        wfe_img = np.zeros((len(x), len(y)))
        aperture1_id = np.zeros((len(x), len(y)), dtype=complex)
        aperture1_ab = np.zeros((len(x), len(y)), dtype=complex)
        aperture1_id_norm = np.zeros((len(x), len(y)))
        aperture1_ab_norm = np.zeros((len(x), len(y)), dtype=complex)
        aperture2 = np.zeros((len(x), len(y)))
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
                
                
                # # specify wafefront error
                wfe_gen = float(wfe(list_wfe, ra, the, D1_2))
                wfe_img[counterx][countery] = float(wfe(list_wfe, ra, the, D1_2))
                
                # define aprture 1
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
                e_field_a2_ab = e_field_a1_ab  * aperture2
                # e_field_a2_plus = e_field_a1_plus * aperture2 #direct
                # e_field_a2_minus = e_field_a1_minus * aperture2 #direct
                # intensity_a2 = abs(e_field_a2.real)**2
                
                # e field in imaging plane
                e_field_id = ifft2(e_field_a2_id)
                e_field_ab = ifft2(e_field_a2_ab)
                e_field_dirty_id = ifft2(e_field_a1_id)
                e_field_dirty_ab = ifft2(e_field_a1_ab)
                
                # e_field_plus = ifft2(e_field_a2_plus) #direct
                # e_field_minus = ifft2(e_field_a2_minus) #direct
                
                # sum and diff of e fields
                e_plus = e_field_id + e_field_ab
                e_minus = e_field_id - e_field_ab
                
                e_plus_dirty = e_field_dirty_id + e_field_dirty_ab
                e_minus_dirty = e_field_dirty_id - e_field_dirty_ab
                
                # calculate intensity in image plane
                intensity_max = abs(e_plus.real)**2
                intensity_min = abs(e_minus.real)**2
                
                
                intensity_max_dirty = abs(e_plus_dirty.real)**2
                intensity_min_dirty = abs(e_minus_dirty.real)**2
                
                # intensity_plus_direct = abs(e_field_plus.real)**2
                # intensity_minus_direct = abs(e_field_minus.real)**2
            
        
        
        # define null
        imax = np.sum(intensity_max)
        imin = np.sum(intensity_min)
        null = imin/imax
        
        # I_init total (|E_1 + E_2|^2)
        iinit_common = np.sum(abs((aperture1_id + aperture1_ab).real)**2)
        
        # define throughput
        throughput = imax/iinit_common
        
        # imaxd = np.sum(intensity_max_dirty)
        # imind = np.sum(intensity_min_dirty)
        # nulld = imind/imaxd
        # print(imax)
        # print(imin)
        # print(imax)
        
        
        
        
        # I_init separat
        # iinit = np.sum(abs(aperture1_id.real)**2) + np.sum(abs(aperture1_ab.real)**2)
        
        # print(np.sum(abs(aperture1_id.real)**2) + np.sum(abs(aperture1_ab.real)**2)
        
        # print(np.sum(abs(aperture1_id.real)**2)) #iinit id
        # print(np.sum(abs(aperture1_ab.real)**2)) # iinit ab
        
        
        
        # print(imax)
        # print(iinit) #iinit sep
        
        print("Common initial intensity: " + str(iinit_common))
        print("Null: " + str(round(null, 10)))
        print("Throughput: " + str(round(throughput, 1)) + " %")
        
        # print(np.sum(intensity_minus_direct)/np.sum(intensity_plus_direct))
        # null = irr_min/irr_max
              
        print("Counters: " + str(counter_pinhole) + " / " + str(counter_rms))
        print("\n")
        
        output_null[counter_pinhole][counter_rms] = null
        output_throughput[counter_pinhole][counter_rms] = throughput
    

'''plotting'''

fig, axs = plt.subplots(1, 2)

# throughput
img1 = axs[0].imshow(output_throughput)
fig.colorbar(img1, ax=axs[0], fraction=0.046, pad=0.04)
axs[0].set_title("Throughout")

# null
img2 = axs[1].imshow(output_null)
fig.colorbar(img2, ax=axs[1], fraction=0.046, pad=0.04)
axs[1].set_title("Null")

# plt.tight_layout(pad=1.5)
# plt.subplots_adjust(wspace=-0.2, hspace=0.71)
plt.savefig("plot_rms.pdf")
plt.show()
