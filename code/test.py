import numpy as np
from scipy.optimize import curve_fit
from peak_fit import *

def voigt_signal(a, p0, g, l, n):
    v0 = voigt(0, g, l)
    return [a * voigt(n - p0, g, l) / v0  for n in range(0, n)]

def poly_signal(args, n):
    return [poly(el, *args) for el in range(n)]

def noise(std, n):
    return list(np.random.randn(n) * std)

n = 5000
a, c, g, l = 1, n / 2,  (n / 20), (n / 30)
yv = voigt_signal(a, c, g, l, n)
yp = poly_signal([1, - 1 / n,  0 / n ** 2, - 0 / n ** 3], n)
std = 0.01
yn = noise(std, n)
y = list(np.sum([yv, yp, yn], axis = 0))
yvp = list(np.sum([yv, yp], axis = 0))

import matplotlib.pyplot as plt
fig = plt.figure(0)
plt.clf()
#plt.plot(np.sum([yv, yn], axis = 0), c = "gray", zorder = 0)
plt.plot(y , c = "gray", zorder = 0)
#yg = peak_detection(y, peak_ratio = 0.1)[0]
x = range(n)
fit, err, snr, fun, yf = fit_voigt(x, y, peak_ratio = 1 / 5, background = 1, log = 1)
plt.plot(yvp, c = "tab:green", zorder = 1)
plt.plot(yf[0], c = "tab:blue", zorder = 1)
plt.xlabel("x")
plt.ylabel("y")
plt.legend(["voigt", "voigt with noise","fitted by peaky"])
fig.set_tight_layout(True) 
plt.show(block = 0)

#print(round(a,1), round(c,1), round(full_width(g, l),1))






