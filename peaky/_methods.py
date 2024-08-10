from scipy.optimize import curve_fit
from scipy.special import wofz as _voigt_scipy
import numpy as np

####### list manipulation #######

reverse = lambda data: data[: : -1]

####### statistical functions #######

def variance(data, order = 0):
    return np.var(data) if order == 0 else np.var([data[i + order] + data[i - order] - 2 * data[i] for i in range(order, len(data) - order)]) / 6

def standard_deviation(data, order = 0):
    return np.sqrt(variance(data, order))

####### curve functions #######

def poly(x, *args):
    pol = 0
    for i in range(len(args)):
        pol += args[i] * np.power(x, i)
    return pol

_kg = 1 / np.sqrt(2 * np.pi)  # 0.398 ...
def gaussian(nu, g):
    if g == 0:
        return 1 * (nu == 0)
    return _kg * g ** -1 * np.exp(- 0.5 * (nu / g) ** 2)

_kl = 1/ np.pi # 0.318 ...
def lorentian(nu, l):
    if l == 0:
        #return 1 if nu == 0 else 0
        return 1 * (nu == 0)
    return _kl * l * (1 / (nu ** 2 + l ** 2))

_kv = 1 / np.sqrt(2) # 0.7071 ...
def voigt(nu, g, l):
    if g == 0:
        return lorentian(nu, l)
    if l == 0:
        return gaussian(nu, g)
    z = _kv * (nu + 1j * l) / g
    w = _voigt_scipy(z)
    return _kg * w.real / g
    
####### averaging  #######

def mean(data):
    return sum(data) / len(data)

def moving_average(data, window = 1, order = 1):
    l = len(data)
    data = continuation(data, window, order)
    g = np.array([1 for i in range(2 * window + 1)])
    return list(np.convolve(data, g, mode = "same")[window: window + l] / (2 * window + 1))

def gaussian_average(data, sigma = 1, n = 3, order = 0):
    l = len(data)
    ns = n * sigma
    data = continuation(data, ns, order)
    g = np.array([gaussian(i, sigma) for i in range(-ns, ns + 1)])
    return list(np.convolve(data, g, mode = "same")[ns: ns + l] / np.sum(g))

####### data manipulation #######

def closest(data, value):
    data = np.array(data) - value
    data = np.where(np.diff(np.sign(data)))[0]
    return list(data)

def derivative(data, order = 1):
    if order == 0:
        return  data
    first = data[1] - data[0]
    last = data[-1] - data[-2]
    d = [first] + [(data[i + 1] - data[i - 1]) / 2 for i in range(1, len(data) - 1)] + [last]
    return derivative(d, order - 1)

def continuation(data, n = 1, order = 1):
    if n == 0:
        return data
    l = len(data)
    n_new = n + order
    x1 = [i for i in range(len(data)) if i in range(0, n_new)]
    x2 = [i for i in range(len(data)) if i in list(range(l - n_new, l))]
    y1 = [data[i] for i in range(len(data)) if i in x1]
    y2 = [data[i] for i in range(len(data)) if i in x2]
    fit1 = reverse(np.polyfit(x1, y1, order))
    fit2 = reverse(np.polyfit(x2, y2, order))
    fit_data1 = [poly(i, *fit1) for i in range(-n, 0)]
    fit_data2 = [poly(i, *fit2) for i in range(l, l + n)]
    data =  fit_data1 + list(data) + fit_data2
    return data

def detrend(data, order = 2, window = 1):
    if window == 0:
        return detrend(data, 1, 1)
    l = len(data)
    pos = [i for i in range(l) if i < window or l - 1 - i < window]
    data_l = [data[i] for i in pos]
    fit_l = reverse(np.polyfit(pos, data_l, order))
    data_l = [poly(i, *fit_l) for i in range(l)]
    mean_l = mean(data_l)
    data = [data[i] - data_l[i] for i in range(l)]
    return data