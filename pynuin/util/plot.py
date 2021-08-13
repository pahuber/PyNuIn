#%%
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft2, fftshift
from pynuin.main.zernike import wfe, noll_index
from pynuin.main.optics import aperture


def plot_zernike_wfe(indices_coeffs, name_plain=None):
    '''
    Quickly plots a sum of given Zernike polynomial terms.

            Parameters:
                    indices_coeffs (list): List containing tuples of the kind (j, coeff), where j is the Zernike poylnomial index according to Noll's convention

            Returns:
                    -
    '''
    
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
                wfe_img[countery][counterx] = float(wfe(indices_coeffs, rho, theta, 1))
            
    # plot image
    plt.imshow(wfe_img)
    plt.colorbar()
    plt.title("Aberrated Wavefront")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xticks([0, len(x)/2, len(x)-1], [-1, 0, 1])
    plt.yticks([0, len(y)/2, len(y)-1], [-1, 0, 1])
        
    if name_plain is not None:
        plt.imsave(name_plain, wfe_img)
        
    plt.show()
    

def plot_zernike_distribution(indices_coeffs, normalized=False, name=None, title=None):
    '''
    Quickly plots the distribution of given Zernike polynomial terms.

            Parameters:
                    indices_coeffs (list): List containing tuples of the kind (j, coeff), where j is the Zernike poylnomial index according to Noll's convention
                    normalized (boolean): If true, coefficient values are normalized to the value of the first coefficient z_1

            Returns:
                    -
    '''
    
    terms = {}
    
    for el in indices_coeffs:
        j = el[0]
        coeff = el[1]
        
        terms[j] = coeff
    
    lists = sorted(terms.items()) 

    x, y = zip(*lists) 
    
    if normalized:
        y /= max(y)
    
    plt.bar(x, y, color="mediumblue")
    if title is None:
        plt.title("Distribution of Zernike Coefficients")
    else:
        plt.title(title)
        
    plt.xlabel("$j$ ($Z_j$)")
    
    if normalized:
        plt.ylabel("Coefficient [$z_{" + str(indices_coeffs[0][0]) + "}$]")
    else:
        plt.ylabel("Coefficient [m]")
        
    plt.grid()
    
    if name is not None:
        plt.savefig(name, bbox_inches = 'tight', pad_inches = 0)
    
    plt.show()
        
    
def plot_aberrated_psf(indices_coeffs, name=None, name_plain=None):
    
    # define aperture
    a1_ab = aperture(D = 5, 
                 lam = 1,
                 a0 = 1,
                 list_wfe = indices_coeffs,
                 n = 100)
    
    # calculate electric field
    e_pinhole_ab = fftshift(fft2(a1_ab))
    
    # calculate intensity, i. e. psf
    intensity = abs(e_pinhole_ab)**2
    
    # plot
    fig, axs = plt.subplots(1, 1)
    img = axs.imshow(intensity)
    # img.set_clim(1e0, np.amax(intensity))
    fig.colorbar(img, ax=axs, fraction=0.046, pad=0.04)   
    
    if name is not None:
        plt.savefig(name, bbox_inches = 'tight', pad_inches = 0)
        
    if name_plain is not None:
        plt.imsave(name_plain, intensity)
        
    plt.show()