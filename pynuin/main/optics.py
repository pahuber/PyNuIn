import numpy as np
from pynuin.main.zernike import wfe
from scipy.optimize import fsolve
        

# method to create aperture matrices
def aperture(D = 1., 
             lam = 1,
             a0 = 1.,
             list_wfe = None,
             xymin = -5,
             xymax = 5,
             steps = 100):
    '''
    Returns an array containing a circular aperture.

            Parameters:
                    D (float): Aperture diameter in meters
                    lam (int): Wavelength in meters
                    a0 (float): Initial amplitude
                    list_wfe(list): List describing the wavefront error
                    xymin (int): Minimum x and y value of the array
                    xymax (int): Maximum x and y value of the array
                    steps (int): Number of steps to take between xymin and xymax

            Returns:
                    aperture (array): Complex array containing a circular aperture of diameter D
    '''
    
    # definitions
    D_2 = D/2

    # define matrix
    x = np.linspace(xymin, xymax, steps)
    y = np.linspace(xymin, xymax, steps)
    aperture = np.zeros((len(x), len(y)), dtype=complex)
    
    # calculate matrix
    for counterx, elx in enumerate(x):
        for countery, ely in enumerate(y):
            
            # perform transformation to polar coordinates
            rho = np.sqrt(elx**2+ely**2)
            theta = np.arctan2(ely, elx)
            
            
            # specify wafefront error
            if list_wfe is not None:
                wfe_gen = float(wfe(list_wfe, rho, theta, D_2))
            else:
                wfe_gen = 0
            
            # define aprture 1
            aperture[counterx][countery] = a0*np.heaviside(D_2-rho, 1)*np.heaviside(rho, 1)*np.exp(-2*np.pi*1j*wfe_gen/lam)
    
    return aperture


# method to normalize two apertures, such that the initial intensity of the sum of the two e-field is intensity_init
def normalize_apertures(a1,
                        a2,
                        intensity_init):

    intensity_init_check = 0
    counter = 0
    a0_imag_new = 0
    
    while round(intensity_init_check, 3) != float(intensity_init):
        
        if counter != 0:
            print("Renormalizing...")
                
        # normalize this amplitude to unit intensity
        aperture_total  = a1 + a2
        amplitude_init = np.sum(aperture_total)
        
        # choose initial amplitude imaginary part because of degeneracy
        if counter == 0:
            a0_imag = amplitude_init.imag
        else:
            a0_imag = a0_imag_new
        
        # define function to find roots for
        def func(a0_real):
            
            func = 0
            
            for el1 in aperture_total:
                for el2 in el1:
                    z_real = el2.real
                    z_imag = el2.imag
                    
                    func += abs(a0_real*z_real - a0_imag*z_imag)**2
                    
            return func - intensity_init
        
        # solve for roots, i. e. real part of amplitude a0
        a0_real = fsolve(func, 1)[0]
        
        def func(a0_imag_new):
            
            func = 0
            
            for el1 in aperture_total:
                for el2 in el1:
                    z_real = el2.real
                    z_imag = el2.imag
                    
                    func += abs(a0_real*z_real - a0_imag_new*z_imag)**2
                    
            return func - intensity_init
        
        # solve for roots, i. e. imaginary part of amplitude a0
        a0_imag_new = fsolve(func, 1)[0]
        
        # define full complex amplitude a0
        a0_total = a0_real + a0_imag_new*1j
        
        # check normalization
        
        intensity_init_check = np.sum(abs(((a0_total*a1+a0_total*a2).real)**2))
        # print(intensity_init_check)
        counter += 1
    
    # print(intensity_init_check)
    
    if counter != 1:
        print("Normalization succesfull")
    
    return a0_total*a1, a0_total*a2, a0_total