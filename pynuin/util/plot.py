#%%
import numpy as np
import matplotlib.pyplot as plt
from pynuin.main.zernike import wfe, noll_index


def plot_zernike_wfe(indices_coeffs):
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
                wfe_img[counterx][countery] = float(wfe(indices_coeffs, rho, theta, 1))
            
    # plot image
    plt.imshow(wfe_img)
    plt.colorbar()
    plt.title("Aberrated Wavefront")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.xticks([0, len(x)/2, len(x)-1], [-1, 0, 1])
    plt.yticks([0, len(y)/2, len(y)-1], [-1, 0, 1])
    plt.show()
    

def plot_zernike_distribution(indices_coeffs):
    '''
    Quickly plots the distribution of given Zernike polynomial terms.

            Parameters:
                    indices_coeffs (list): List containing tuples of the kind (j, coeff), where j is the Zernike poylnomial index according to Noll's convention

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
    
    plt.bar(x, y)
    plt.title("Distribution of Zernike Polynomial Terms")
    plt.xlabel("Zernike Term")
    plt.ylabel("Coefficient [m]")
    plt.grid()
    plt.show()
        