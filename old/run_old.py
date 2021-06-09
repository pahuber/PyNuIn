#%%

import numpy as np
import matplotlib.pyplot as plt

from pynuin.main.model import wfe

from scipy.special import j1
from numpy.fft import fft2, fftshift, ifft2

from matplotlib.colors import LogNorm



# prepare matrices
x = np.arange(-10, 10, 0.2)
y = np.arange(-10, 10, 0.2)

test = np.zeros((len(x), len(y)))

img = np.zeros((len(x), len(y)))
amp_plus = np.zeros((len(x), len(y)), dtype=complex)
amp_minus = np.zeros((len(x), len(y)), dtype=complex)
i_max = np.zeros((len(x), len(y)))
i_min = np.zeros((len(x), len(y)))
null = np.zeros((len(x), len(y)))




# # test e field
# list_wfe = [(1, 1, 250e-9), (1, -1, 250e-9)]
# wfe_gen = wfe(list_wfe)

# aperture_fct = circular_aperture()

# print(aperture_fct)

# print(wfe_gen)


# e_field = e_id(e0=1e-8,
#           t=1,
#           x=1e-2,
#           y=1e-2,
#           z=1e-1,
#           D=1e-2, 
#           lam=450e-9, 
#           rho_max=1e-1, #2*D
#           c=3e8)


# e_field = e_id_polar(t=1, 
#          q=1e-2, 
#          z=1e-1,
#          e0=1e-8,
#          D=1e-2, 
#          lam=450e-9,
#          c=3e8)

# print(e_field)


#%%
# draw e field

# for counterx, elx in enumerate(x):
#     for countery, ely in enumerate(y):
        
#         # print(counterx)
#         # print(countery)
        
#         # perform transformation to polar coordinates
#         ra = np.sqrt(elx**2+ely**2)
#         the = np.arctan2(ely, elx)

#         # if ra < 1:
#         test[counterx][countery] = e_id_polar(t=0, q=ra, z=1e6, e0=1e-8, D=1e-3, lam=1.1e-5, c=3e8)
        
#         # test[counterx][countery] = j1(ra/np.sqrt(ra**2+1))
#         # e_id(1e-8,t=1,x=elx,y=ely,z=1, D=1e-2, lam=450e-9, rho_max=1e-1, c=3e8)
#         print(test[counterx][countery])
            
        


# # plot matrices        
# plt.imshow(test)
# # plt.title("Null Depth")
# # plt.clim(-1e-18, 0)
# plt.colorbar()
# # plt.xticks([0, 20, 40], [-1, 0, 1])
# # plt.yticks([0, 20, 40], [-1, 0, 1])
# plt.show()


#%%

# draw aperture
for counterx, elx in enumerate(x):
    for countery, ely in enumerate(y):
        
        # print(counterx)
        # print(countery)
        
        # perform transformation to polar coordinates
        ra = np.sqrt(elx**2+ely**2)
        the = np.arctan2(ely, elx)
        
        a0 = 1
        D_2 =0.4
        lam = 1e-5
        
        if 4 < elx and elx < 6:
            if 4 < ely and ely < 6:
                test[counterx][countery] = 1
            
        
        
        # if ra < D_2:
            # test[counterx][countery] = 1
        
        # list_wfe = [(1, 1, 0.5e-5), (1, -1, 0.5e-5)]
        # # wfe_gen = float(wfe_cart(list_wfe, elx, ely, max(x), max(y)))
        # wfe_gen = float(wfe(list_wfe, ra, the))

        # # ideal e field
        # # test[counterx][countery] = a0*np.heaviside(D_2-ra, 1)*np.heaviside(ra, 1)
        
        # # aberrated e field
        # test[counterx][countery] = a0*np.heaviside(D_2-ra, 1)*np.heaviside(ra, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)
        
        # # sum
        # amp_plus[counterx][countery] = a0*np.heaviside(D_2-ra, 1)*np.heaviside(ra, 1) + a0*np.heaviside(D_2-ra, 1)*np.heaviside(ra, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)
            
        # # difference
        # amp_minus[counterx][countery] = a0*np.heaviside(D_2-ra, 1)*np.heaviside(ra, 1) - a0*np.heaviside(D_2-ra, 1)*np.heaviside(ra, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)


# test2 = abs(fftshift(fft2(amp_plus)).real)**2
# test3 = abs(fftshift(fft2(amp_minus)).real)**2

# null = test3/test2

# for el1 in test2:
#     for el2 in el1:
#         print(el2.real)

# print(test2.shape)        


# plot matrices


plt.imshow(fftshift(abs(ifft2(test).real))**2)
# print(null)
# plt.title("Null Depth")
# plt.clim(0, 10)
plt.colorbar()
# plt.xticks([0, 20, 40], [-1, 0, 1])
# plt.yticks([0, 20, 40], [-1, 0, 1])
plt.show()

#%%

# draw zernike

# # calculate matrices
# for counterx, elx in enumerate(x):
#     for countery, ely in enumerate(y):
                
#         # perform transformation to polar coordinates
#         ra = np.sqrt(elx**2+ely**2)
#         the = np.arctan2(ely, elx)
#         print(the)

#         if ra < 1:
#             # define wavefront error        
#             wfe1 = float(wfe_basic(piston=0,
#                              tiltx=250e-9,
#                              tilty=250e-9,
#                              defocus=0,
#                              rho=ra,
#                              theta=the))
            
#             # define general wavefront error
#             list_wfe = [(1, 1, 250e-9), (1, -1, 250e-9)]
#             wfe_gen = float(wfe(list_wfe, rho=ra, theta=the))
            
#             # define amplitudes
#             a_id = amp_id(1, t=0, lam=450e-9)
#             a_ab = amp_ab(1, wfe_gen, t=0, lam=450e-9, rho=ra, theta=the)
            
#             # define matrices
#             amp_plus[counterx][countery] = a_id + a_ab
#             amp_minus[counterx][countery] = a_id - a_ab 
#             i_max[counterx][countery] = abs(a_id + a_ab)**2
#             i_min[counterx][countery] = abs(a_id - a_ab)**2
#             null[counterx][countery] = (abs(a_id - a_ab)**2)/(abs(a_id + a_ab)**2)
            
#             img[counterx][countery] = z(1, 1, ra, the) + z(1, -1, ra, the)
        
        
        


# # plot matrices        
# plt.imshow(img)
# # plt.title("Null Depth")
# plt.colorbar()
# # plt.xticks([0, 20, 40], [-1, 0, 1])
# # plt.yticks([0, 20, 40], [-1, 0, 1])
# plt.show()


# test