import numpy as np
from pynuin.main.zernike import wfe
from scipy.optimize import fsolve
from mpmath import findroot
        

def aperture(D = 1., 
             lam = 1,
             a0 = 1.,
             list_wfe = None,
             n = 100):
    '''
    Returns a numpy array containing a circular aperture, possibly including wavefront aberrations.

            Parameters:
                    D (float): Aperture diameter in meters
                    lam (int): Wavelength in meters
                    a0 (float): Initial amplitude
                    list_wfe(list): List describing the wavefront error
                    n (int): Number to specify shape of array, i. e. n x n

            Returns:
                    (array): Complex array containing a circular aperture of diameter D
    '''
    
    # create 1-dimensional x and y arrays
    i = np.arange(1, n)
    x = i - n / 2
    y = n / 2 - i
        
    # convert to 2-dimensional arrays
    Y = y[:, np.newaxis]
    X = x[np.newaxis, :]
    
    # define radius
    R = D/2
    
    # define ideal aperture
    aperture = a0*(X**2 + Y**2 < R**2).astype(complex)
    
    # add wavefront error if a list of terms is given
    if list_wfe is not None:
        for counterx, elx in enumerate(x):
            for countery, ely in enumerate(y):
                aperture[countery, counterx] *= np.exp(-2*np.pi*1j*float(wfe(list_wfe, np.sqrt(elx**2+ely**2), np.arctan2(ely, elx), R))/lam)
        
    return aperture