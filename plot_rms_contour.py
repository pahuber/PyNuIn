#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe, get_coeff_from_rms, get_distribution_from_rms
from pynuin.main.optics import aperture
from pynuin.util.plot import plot_zernike_distribution
from numpy.fft import fft2, ifft2, fftshift
from matplotlib.colors import LogNorm
from matplotlib import ticker
from scipy.optimize import fsolve


'''computational specifications'''
n = 100
size = n**2


'''physical specifications'''
lam = 1 #1e-5 #m
D1 = 1*20 #0.01 #m
# list_wfe = [(2, 0), (2, 2), (2, -2), (3, 1), (3, -1), (3, 3), (3, -3), (4, 0), (4, 2), (4, -2), (4, 4), (4, -4), (5, 5), (5, 3), (5, 1), (5, -5), (5, -3), (5, -1)]
list_wfe = []
for i in range(5, 5+1, 1):
    list_wfe.append((i,))

pinholes = np.arange(0.2*lam/D1*n, 5.2*lam/D1*n, 0.2*lam/D1*n)
rmss = np.arange(0.001*lam, 0.026*lam, 0.001*lam)

print(len(pinholes))
print(len(rmss))
print(len(pinholes) * len(rmss))

# define output arrays
output_null = np.zeros((len(rmss), len(pinholes)))
output_throughput = np.zeros((len(rmss), len(pinholes)))


'''start loop'''
for counter_rms, rms in enumerate(rmss):
    
    '''rms to zernike coefficient calculations'''
    list_wfe2 = list_wfe.copy()
    
    # equally distributed
    coeff = get_coeff_from_rms(rms, list_wfe2)
    for index in range(len(list_wfe)):
        list_wfe2[index] += (coeff,)
    
    # unequally distibuted, decreasing values
    # list_wfe2 = get_distribution_from_rms(rms, list_wfe2, rms/40) # 40, such that the slope is that no coefficients are negative
    # 
    # plot distributions
    # plot_zernike_distribution(list_wfe2, normalized=True)
    
    
    '''create apertures'''
    # create ideal aperture
    a1_id = aperture(D = D1, 
                 lam = lam,
                 a0 = 1,
                 list_wfe = None,
                 n = n)
    
    # create aberrated aperture
    a1_ab = aperture(D = D1, 
                 lam = lam,
                 a0 = 1,
                 list_wfe = list_wfe2,
                 n = n)
    
    
    # normalize apertures
    # print("Normalizing")
    # a1_id, a1_ab, _ = normalize_apertures2(a1_id, a1_ab, 1)
    # print("Ended normalizing")
    
    
    for counter_pinhole, pinhole in enumerate(pinholes):
            
        # create aperture in pinhole plane
        a2 = aperture(D = pinhole, 
              lam = lam,
              a0 = 1,
              list_wfe = None,
              n = n)
        
        
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
        intensity_max = abs(e_plus)**2
        intensity_min = abs(e_minus)**2
        
        # define null
        imax = np.sum(intensity_max)
        imin = np.sum(intensity_min)
        null = imin/imax
        # print(null)
        
        # calculate initial common intensity, should equal intensity_init from above, i. e. I_init total (|E_1 + E_2|^2)
        iinit_common = np.sum(abs((a1_id + a1_ab))**2)
        
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
tcontours = [0.5, 0.7, 0.80, 0.85, 0.90]
tcontours_full = [0.2, 0.5, 0.7, 0.8, 0.85, 0.90, 1.0]
# tcontours = [0.1, 0.15, 0.2]
ncontours = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
cmap = "viridis"
cmap_r = "viridis_r"
tcolor="white"
ncolor="black"

# throughput
fig1, axs1 = plt.subplots(1, 1)
img1 = axs1.imshow(output_throughput, extent=extent, aspect=aspect, origin="lower", cmap=cmap)


cs1 = axs1.contour(output_throughput, tcontours_full, colors=tcolor, extent=extent, linewidths=1.5)
axs1.clabel(cs1, cs1.levels, inline=True, fmt="%1.2f", fontsize=10)

cs12 = axs1.contour(output_null, ncontours, colors=ncolor, extent=extent, linewidths=1.5)
axs1.clabel(cs12, cs12.levels, inline=True, fmt="%.0E", fontsize=10)

fig1.colorbar(img1, ax=axs1, fraction=0.046, pad=0.04)
axs1.set_xticks(pinholes)
axs1.set_xticklabels((pinholes/lam/n*D1).round(2))
axs1.set_yticks(rmss)
axs1.set_title("Throughput (D2)")
plt.xlabel("Pinhole Diameter [$\lambda/D$]")
plt.ylabel("RMS [$\lambda$]")
plt.xticks(rotation=90, fontsize=7)
plt.yticks(fontsize=7)

fig1.savefig("output/rms_throughput.pdf", bbox_inches = 'tight', pad_inches = 0)

# null
fig2, axs2 = plt.subplots(1, 1)
img2 = axs2.imshow(output_null, norm=LogNorm(), extent=extent, aspect=aspect, origin="lower", cmap=cmap_r)

cs2 = axs2.contour(output_throughput, tcontours_full, colors=tcolor, extent=extent, linewidths=1.5)
axs2.clabel(cs2, cs2.levels, inline=True, fmt="%1.2f", fontsize=10)

cs22 = axs2.contour(output_null, ncontours, colors=ncolor, extent=extent, linewidths=1.5)
axs2.clabel(cs22, cs22.levels, inline=True, fmt="%.0E", fontsize=10)

# axs2.fill_between(pinholes, cs2, cs22)

fig2.colorbar(img2, ax=axs2, fraction=0.046, pad=0.04)
axs2.set_xticks(pinholes)
axs2.set_xticklabels((pinholes/lam/n*D1).round(2))
axs2.set_yticks(rmss)
axs2.set_title("Null Depth (D2)")
plt.xlabel("Pinhole Diameter [$\lambda/D$]")
plt.ylabel("RMS [$\lambda$]")
plt.xticks(rotation=90, fontsize=7)
plt.yticks(fontsize=7)

fig2.savefig("output/rms_null.pdf", bbox_inches = 'tight', pad_inches = 0)

plt.show()

# plot distirbution
plot_zernike_distribution(list_wfe2, normalized=True, name="output/rms_distro.pdf", title="Unequally Distributed Coefficients (D2)")




# plot contour


fig, ax = plt.subplots()

# throughput axis
ax2 = ax.twinx()

# throughput contour for axis scaling -------

tcontours_mapto_pinhole = []

for counter, contour in enumerate(cs1.collections):

    X = []
    Y = []
    label = tcontours_full[counter]
    # print(label)

    
    for path in contour.get_paths():
        v = path.vertices
        x = v[:,0]
        y = v[:,1]
        
        X += x.tolist() #pinhole
        Y += y.tolist()
        
    Xmean = np.mean(X)
    tcontours_mapto_pinhole.append(Xmean)
    
# ---------
    

# # create equally space throughput axis by linearly interpolating values
# throughput_interpolation = []

# for counter in range(len(tcontours_mapto_pinhole)-2):
#     if counter != len(tcontours_mapto_pinhole)-1:
#         plow = tcontours_mapto_pinhole[counter]
#         phigh = tcontours_mapto_pinhole[counter+1]
#         tlow = tcontours_full[counter]
#         thigh = tcontours_full[counter+1]
        
#         pnum = pinholes.tolist().index(phigh.round(0)) - pinholes.tolist().index(plow.round(0))
        
#         list_temp = np.linspace(tlow, thigh + (thigh-tlow)/pnum, pnum)
        
        
#         throughput_interpolation += list_temp.tolist()
    
    
# # -------

# print(len(pinholes))
# print(len(throughput_interpolation))
# print(throughput_interpolation)



# null depth contours for plotting
for counter, contour in enumerate(cs12.collections):

    X = []
    Y = []
    label = "N = " + str("{:.0e}".format(ncontours[counter]))

    
    for path in contour.get_paths():
        v = path.vertices
        x = v[:,0]
        y = v[:,1]
        
        X += x.tolist() #pinholes ~~
        Y += y.tolist()

    ax.plot(Y, X, label=label)
    ax2.plot(Y, X)
        
plt.title("Pinhole Diameter/Throughput vs. RMS (D2)")
ax.set_xlabel("RMS [$\lambda$]")
ax.set_ylabel("Pinhole Diameter [$\lambda/D$]")
ax.set_yticks(pinholes)
ax.set_yticklabels((pinholes/lam/n*D1).round(2), fontsize=7)
ax.legend()
ax.grid(axis="x")

ax2.set_ylabel("Throughput")
ax2.set_yticks(tcontours_mapto_pinhole)
ax2.set_yticklabels(tcontours_full, fontsize=9)
ax2.grid()

plt.savefig("output/rms_contour.pdf", bbox_inches = 'tight', pad_inches = 0)