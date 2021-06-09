#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe, get_coeff_from_rms
from pynuin.main.optics import aperture, normalize_apertures
from numpy.fft import fft2, ifft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker
from scipy.optimize import fsolve


'''definitions'''
# define input
lam = 1 #1e-5 #m
D1 = 1 #0.01 #m
list_wfe = [(2, 0), (2, 2), (2, -2), (3, 1), (3, -1), (3, 3), (3, -3), (4, 0), (4, 2), (4, -2), (4, 4), (4, -4), (5, 5), (5, 3), (5, 1), (5, -5), (5, -3), (5, -1)]
list_wfe = [(2, 2)]
pinholes = np.arange(1*lam/D1, 3*lam/D1, 0.2*lam/D1)
rmss = np.arange(0.005*lam, 0.055*lam, 0.005*lam)
xymin = -5
xymax = 5
steps = 500

# print(len(pinholes))
# print(len(rmss))

print(len(pinholes) * len(rmss))

# define output arrays
output_null = np.zeros((len(rmss), len(pinholes)))
output_throughput = np.zeros((len(rmss), len(pinholes)))

for counter_rms, rms in enumerate(rmss):
    
    
    '''rms to zernike coefficient calculations'''
    list_wfe2 = list_wfe.copy()
    coeff = get_coeff_from_rms(rms, list_wfe2)
    for index in range(len(list_wfe)):
        list_wfe2[index] += (coeff,)
    
    
    '''create apertures'''
    # create ideal aperture
    a1_id = aperture(D = D1, 
                  lam = lam,
                  a0 = 1,
                  list_wfe = None,
                  xymin = xymin,
                  xymax = xymax,
                  steps = steps)
    
    # create aberrated aperture
    a1_ab = aperture(D = D1, 
                  lam = lam,
                  a0 = 1,
                  list_wfe = list_wfe2,
                  xymin = xymin,
                  xymax = xymax,
                  steps = steps)
    
    
    # normalize apertures
    a1_id, a1_ab, _ = normalize_apertures(a1_id, a1_ab, 1)
    
    
    for counter_pinhole, pinhole in enumerate(pinholes):
            
        # create aperture in pinhole plane
        a2 = aperture(D = pinhole, 
                      lam = lam,
                      a0 = 1,
                      list_wfe = None,
                      xymin = xymin,
                      xymax = xymax,
                      steps = steps)
        
        
        '''calculations'''
        # e field in pinhole plane
        e_pinhole_id = fftshift(fft2(a1_id))
        e_pinhole_ab = fftshift(fft2(a1_ab))
        
        # e field in final plane
        e_pinhole_filtered_id = e_pinhole_id * a2
        e_pinhole_filtered_ab = e_pinhole_ab  * a2
        e_final_id = ifft2(e_pinhole_filtered_id)
        e_final_ab = ifft2(e_pinhole_filtered_ab)
        
        # sum and diff of e fields
        e_plus = e_final_id + e_final_ab
        e_minus = e_final_id - e_final_ab
        
        # calculate intensity in final plane
        intensity_max = abs(e_plus.real)**2
        intensity_min = abs(e_minus.real)**2
        
        # define null
        imax = np.sum(intensity_max)
        imin = np.sum(intensity_min)
        null = imin/imax
        # print(null)
        
        # calculate initial common intensity, should equal intensity_init from above, i. e. I_init total (|E_1 + E_2|^2)
        iinit_common = np.sum(abs((a1_id + a1_ab).real)**2)
        
        # define throughput
        throughput = imax/iinit_common
        
        # print("Common initial intensity: " + str(iinit_common))
        # print("Null: " + str(round(null, 10)))
        # print("Throughput: " + str(round(imax/iinit_common * 100, 1)) + " %")
        
        output_null[counter_rms][counter_pinhole] = null
        output_throughput[counter_rms][counter_pinhole] = throughput
        
        # if counter_rms == 2 and counter_pinhole == 0:
            # print(rms)
            # print(pinhole)
        
# output_null[2][0] = 10000
    

'''plotting'''
pmin = np.min(pinholes)
pmax = np.max(pinholes)
rmin = np.min(rmss)
rmax = np.max(rmss)
ppixel = (pmax-pmin)/len(pinholes)/2
rpixel = (rmax-rmin)/len(rmss)/2
extent = [pmin-ppixel, pmax+ppixel, rmin-rpixel, rmax+rpixel]
aspect = ((pmax-pmin)/output_null.shape[1]) / ((rmax-rmin)/output_null.shape[0])
tcontours = [0.80, 0.85, 0.90]
# tcontours = [0.1, 0.15, 0.2]
ncontours = [1e-6, 1e-5, 1e-4]
cmap = "viridis"
cmap_r = "viridis_r"
tcolor="white"
ncolor="black"

# throughput
fig1, axs1 = plt.subplots(1, 1)
img1 = axs1.imshow(output_throughput, extent=extent, aspect=aspect, origin="lower", cmap=cmap)
cs1 = axs1.contour(output_throughput, tcontours, colors=tcolor, extent=extent, linewidths=1.5)
axs1.clabel(cs1, cs1.levels, inline=True, fmt="%1.2f", fontsize=10)

cs12 = axs1.contour(output_null, ncontours, colors=ncolor, extent=extent, linewidths=1.5)
axs1.clabel(cs12, cs12.levels, inline=True, fmt="%.0E", fontsize=10)

fig1.colorbar(img1, ax=axs1, fraction=0.046, pad=0.04)
axs1.set_xticks(pinholes)
axs1.set_yticks(rmss)
axs1.set_title("Throughput")
plt.xlabel("Pinhole Diameter [$\lambda/D$]")
plt.ylabel("RMS [$\lambda$]")
plt.xticks(rotation=45)

fig1.savefig("rms_throughput.pdf")

# null
fig2, axs2 = plt.subplots(1, 1)
img2 = axs2.imshow(output_null, norm=LogNorm(), extent=extent, aspect=aspect, origin="lower", cmap=cmap_r)

cs2 = axs2.contour(output_throughput, tcontours, colors=tcolor, extent=extent, linewidths=1.5)
axs2.clabel(cs2, cs2.levels, inline=True, fmt="%1.2f", fontsize=10)

cs22 = axs2.contour(output_null, ncontours, colors=ncolor, extent=extent, linewidths=1.5)
axs2.clabel(cs22, cs22.levels, inline=True, fmt="%.0E", fontsize=10)

fig2.colorbar(img2, ax=axs2, fraction=0.046, pad=0.04)
axs2.set_xticks(pinholes)
axs2.set_yticks(rmss)
axs2.set_title("Null Depth")
plt.xlabel("Pinhole Diameter [$\lambda/D$]")
plt.ylabel("RMS [$\lambda$]")
plt.xticks(rotation=45)

fig2.savefig("rms_null.pdf")

plt.show()