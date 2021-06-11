# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 14:05:48 2021

@author: philh
"""


import numpy as np
import pylab as py
from scipy import misc, fftpack

n = 100
I = np.arange(1, n)
x = I - n / 2
y = n / 2 - I


R = 10

X = x[:, np.newaxis]
Y = y[np.newaxis, :]



M = (X**2 + Y**2 < R**2).astype(int)

print(M)

D1 = fftpack.fft2(M)
D2 = fftpack.fftshift(D1)

abs_image = np.abs(D2)


print(np.sum(abs_image**2)/abs_image.size)

py.imshow(abs_image)
py.colorbar()
py.show()