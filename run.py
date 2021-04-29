import numpy as np
import matplotlib.pyplot as plt

from main.model import WFE, A_id, A_ab



# prepare matrices
x = np.arange(-1, 1, 0.05)
y = np.arange(-1, 1, 0.05)
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
        the = abs(np.arctan(ely/elx))
        
        wfe1 = float(WFE(piston=1e-8, tiltx=1e-8, tilty=1e-8, defocus=1e-9, rho=ra, theta=the))
        a_id = A_id(1, t=0, lam=450e-9)
        a_ab = A_ab(1, wfe1, t=0, lam=450e-9, rho=ra, theta=the)
        
        # print("-------------------")
        # print(a_ab)
        # print("-------------------")
        
        
        
        amp_plus[counterx][countery] = a_id + a_ab
        amp_minus[counterx][countery] = a_id - a_ab 
        i_max[counterx][countery] = abs(a_id + a_ab)**2
        i_min[counterx][countery] = abs(a_id - a_ab)**2
        null[counterx][countery] = (abs(a_id - a_ab)**2)/(abs(a_id + a_ab)**2)
        
        # img[counterx][countery] = z(3, -1, ra, the)


# plot matrices        
plt.imshow(null)
plt.colorbar()
plt.show()