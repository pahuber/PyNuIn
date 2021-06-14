#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe


# method to quickly plot an aberrated wavefront
def plot_wfe(list_wfe):
    
    x = np.arange(-1, 1, 0.05)
    y = np.arange(-1, 1, 0.05)
    wfe_img = np.zeros((len(x), len(y)))
    
    for counterx, elx in enumerate(x):
        for countery, ely in enumerate(y):
            
            # perform transformation to polar coordinates
            rho = np.sqrt(elx**2+ely**2)
            theta = np.arctan2(ely, elx)
            
            # calculate image
            if rho <= 1:
                wfe_img[counterx][countery] = float(wfe(list_wfe, rho, theta, 1))
            
    # plot image
    plt.imshow(wfe_img)
    plt.colorbar()
    plt.title("Aberrated Wavefront")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xticks([0, len(x)/2, len(x)-1], [-1, 0, 1])
    plt.yticks([0, len(y)/2, len(y)-1], [-1, 0, 1])
    plt.show()
    
    
# # method to quickly plot an aberrated wavefront
# def plot_wfe_from_rms(rms, list_wfe):
    
#     x = np.arange(-1, 1, 0.5)
#     y = np.arange(-1, 1, 0.5)
#     wfe_img = np.zeros((len(x), len(y)))
    
#     for counterx, elx in enumerate(x):
#         for countery, ely in enumerate(y):
            
#             # perform transformation to polar coordinates
#             rho = np.sqrt(elx**2+ely**2)
#             theta = np.arctan2(ely, elx)
            
#             # calculate image
#             if rho <= 1:
#                 wfe_img[counterx][countery] = float(get_wfe_from_rms(rms, list_wfe, rho, theta))
            
#     # plot image
#     plt.imshow(wfe_img)
#     plt.colorbar()
#     plt.title("Aberrated Wavefront")
#     plt.xlabel("x")
#     plt.ylabel("y")
#     plt.xticks([0, len(x)/2, len(x)-1], [-1, 0, 1])
#     plt.yticks([0, len(y)/2, len(y)-1], [-1, 0, 1])
#     plt.show()