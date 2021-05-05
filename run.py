import numpy as np
import matplotlib.pyplot as plt

from pynuin.main.model import wfe, wfe_basic, amp_id, amp_ab, z



# prepare matrices
x = np.arange(-1, 1.05, 0.05)
y = np.arange(-1, 1.05, 0.05)
img = np.zeros((len(x), len(y)))
amp_plus = np.zeros((len(x), len(y)))
amp_minus = np.zeros((len(x), len(y)))
i_max = np.zeros((len(x), len(y)))
i_min = np.zeros((len(x), len(y)))
null = np.zeros((len(x), len(y)))



# calculate matrices
for counterx, elx in enumerate(x):
    for countery, ely in enumerate(y):
                
        # perform transformation to polar coordinates
        ra = np.sqrt(elx**2+ely**2)
        the = np.arctan2(ely, elx)
        print(the)

        # define wavefront error        
        wfe1 = float(wfe_basic(piston=0,
                         tiltx=250e-9,
                         tilty=250e-9,
                         defocus=0,
                         rho=ra,
                         theta=the))
        
        list_wfe = [(1, 1, 250e-9), (1, -1, 250e-9)]
        wfe_gen = float(wfe(list_wfe, rho=ra, theta=the))
        
        # define amplitudes
        a_id = amp_id(1, t=0, lam=450e-9)
        a_ab = amp_ab(1, wfe_gen, t=0, lam=450e-9, rho=ra, theta=the)
        
        # define matrices
        amp_plus[counterx][countery] = a_id + a_ab
        amp_minus[counterx][countery] = a_id - a_ab 
        i_max[counterx][countery] = abs(a_id + a_ab)**2
        i_min[counterx][countery] = abs(a_id - a_ab)**2
        null[counterx][countery] = (abs(a_id - a_ab)**2)/(abs(a_id + a_ab)**2)
        
        img[counterx][countery] = z(1, 1, ra, the) + z(1, -1, ra, the)


# plot matrices        
plt.imshow(img)
# plt.title("Null Depth")
plt.colorbar()
plt.show()


# test