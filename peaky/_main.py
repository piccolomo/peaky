from ._methods import *
from ._string import *
from ._doc import *


####### peak detection  #######

_c2 = 2 * np.sqrt(2 * np.log(2))
_k1 = 0.5346
_k2 = 0.2166
def full_width(g, l): #fwhm
    g = _c2 * g
    l = 2 * l
    width = _k1 * l + np.sqrt(_k2 * l ** 2 + g ** 2)
    return width
    
def full_width_error(g, l, dg, dl):
    g = _c2 * g
    l = 2 * l
    dg = _c2 * dg
    dl = 2 * dl
    error_g = dg * g / np.sqrt(g ** 2 + _k2 * l ** 2)
    error_l = dl * (_k1 + _k2 * l / np.sqrt(g ** 2 + _k2 * l ** 2))
    return np.sqrt(error_l ** 2 + error_g ** 2)
    
def peak_guess(data, log = True):
    m = max(data)
    p_center = round(int(mean(closest(data, m))))
    a = data[p_center]
    p_left = closest(data[0 : p_center], m / 2)
    p_right = closest(data[p_center : ], m / 2)
    if p_left != [] and p_right != []:
        p_left = p_left[-1]
        p_right = p_center + p_right[-1]
        return a, p_center, (p_right - p_left)
    else:
        if log:
            print(_set_color("warning: initial guess not trusted", "orange"))
        return None, None, None

def peak_detection(data, peak_ratio = 0.1, log = True):
    l = len(data)
    window = int(round(0.1 * l))
    data_d = detrend(data, 3, window)
    sigma = int(round(peak_ratio * l / 8))
    data_g = gaussian_average(data_d, sigma, 2)
    
    a, c, width = peak_guess(data_g, log)
    if a != None and width != None:
        if a < 2 * standard_deviation(data, 1) or width > int(round(0.9 * l)):
            if log:
                print(set_color("warning: initial guess not trusted", "orange"))
            a, c, width = None, None, None
    if width != None:
        kg = voigt(0, width/(2 * _c2), width / (2 * 2))
        if width > _c2 * sigma:
            width = np.sqrt(width ** 2 - (_c2 * sigma) ** 2)
        k = voigt(0, width/(2 * _c2), width / (2 * 2))
        a = a * k / kg
    return a, c, width

def fit_voigt(x, y, peak_ratio = 0.1, background = 1, a = None, p = None, g = None, l = None, log = True):   
    x, y = np.array(x), np.array(y)
    
    std = standard_deviation(y, 1)
    a_min, a_max = 0, 2 * (max(y) - min(y))
    xc_min, xc_max = min(x), max(x)
    g_min, g_max = 0, 1 * (max(x) - min(x))
    l_min, l_max = g_min, g_max

    ai, ci, width = peak_detection(y, peak_ratio, log)
    n = len(y)
    if ai == None:
        ai, ci, width = std, int(round(n / 2)), int(round(0.9 * n))
    xc = x[ci]
    width = (max(x) - min(x)) * width / n
    width_g = width / 2
    width_l = width / 2
    gi, li =  width_g / _c2, width_l / 2

    fit_l = []
    if background != None and background >= 0:
        left , right = ci - int(round(width / 2)), ci + int(round(width / 2))
        x_l = [x[i] for i in range(n) if  i < left or i > right]
        y_l = [y[i] for i in range(n) if  i < left or i > right]
        fit_l = reverse(np.polyfit(x_l, y_l, background))
        fit_l = [fit_l[i] * n ** i for i in range(background + 1)]

    voigt_model = lambda nu, a, nu_0, g, l: a * voigt(nu - nu_0, g, l) / voigt(0, g, l)
    background_model = lambda nu, *args: poly(nu / n, *args)
    signal_model = lambda nu, a, nu_0, g, l, *args: voigt_model(nu, a, nu_0, g, l) + background_model(nu, *args)
    
    iniz = [ai, xc, gi, li] + fit_l
    lb = [a_min, xc_min, g_min , l_min]  # lower bounds
    ub = [a_max, xc_max, g_max, l_max]  # upper bounds
    lb += [-np.inf] * (len(iniz) - 4)
    ub += [np.inf] * (len(iniz) - 4)

    var = np.array([a, p, g, l] + [None] * (len(iniz) - 4))
    iniz = np.array(iniz)[var == None]
    lb = np.array(lb)[var == None]
    ub = np.array(ub)[var == None]
    bounds = (tuple(lb), tuple(ub))

    fit_model = lambda nu, *args: signal_model(nu, *_old_var(args, var))
    
    try:
        fit, cov = curve_fit(fit_model, x, y, p0 = iniz, maxfev = 5 * 10 ** 4, method = "trf", bounds = bounds)
        err = np.sqrt(np.diag(cov))
    except:
        if log:
            print(_set_color("warning: fit failed", "red"))
        fit = iniz
        err = np.array(ub) - np.array(lb)
    fit, err = list(fit), list(err)
    
    iniz = _old_var(iniz, var)
    fit = _old_var(fit, var)
    err = _old_var(err, [el if el == None else 0 for el in var])
    
    functions = [lambda nu: signal_model(nu, *fit), lambda nu: voigt_model(nu, *fit[:4]), lambda nu: background_model(nu, *fit[4:])]

    data = [f(x) for f in functions]
    
    iniz2, fit2, err2 = iniz[:], fit[:], err[:]
    iniz2.insert(4, full_width(iniz[2], iniz[3]))
    fit2.insert(4, full_width(fit[2], fit[3]))
    err2.insert(4, full_width_error(fit[2], fit[3], err[2], err[3]))

    err_f = np.sqrt(mean([(y[i] - data[0][i]) ** 2 for i in range(n)]))
    snr = fit[0] / err_f

    if log:
        print_fit(fit2, err2, iniz2)
        print_snr(snr)
        
    return fit2, err2, snr, functions, data

def _old_var(new_var, var): #where var is None, old_
    old = [0] * len(var)
    j = 0
    for i in range(len(var)):
        if var[i] == None:
            old[i] = new_var[j]
            j += 1
        else:
            old[i] = var[i]
    return old
    
fit_gaussian = lambda x, y, peak_ratio = 0.1, background = 1, a = None, p = None, g = None, log = True: fit_voigt(x, y, peak_ratio, background, a = a, p = p, g = g, l = 0, log = log)

fit_lorentian = lambda x, y, peak_ratio = 0.1, background = 1, a = None, p = None, l = None, log = True: fit_voigt(x, y, peak_ratio, background, a = a, p = p, g = 0, l = l, log = log)

def print_fit(fit, err, iniz):
    l = len(iniz)
    header = ["amplitude", "position", "gaussian", "lorentian", "fwhm"]
    header += ["background" + str(i) for i in range(0, l - 5)]
    header = list_to_string(header, "right")

    per = []
    for i in range(len(header)):
        if err[i] != 0 and fit[i] == 0:
            per.append("inf")
        elif err[i] == 0 and fit[i] == 0:
            per.append("----")
        else:
            per_el = 1000 * abs(err[i] / fit[i])
            if per_el <= 100:
                per.append(per_el)
            else:
                per.append(">100")
            
    err = [err[i] if abs(err[i]) <= abs(fit[i]) else "----" for i in range(len(err))]
    fit = list_to_string(fit, "left", dec = 1)
    err = list_to_string(err, "left", dec = 2)
    per = list_to_string(per, "left", dec = 1)
    iniz = list_to_string(iniz, "left", dec = 1)

    for i in range(len(header)):
        print(header[i] + " " + fit[i] + " ± " + err[i] + " (" + per[i] + " ‰) [initial guess: " + iniz[i] + "]")
        if i == 3:
            pass
            #print()
    #print()

def run_test():
    n = 1000
    std = 0.01
    a, p0, g, l = 1, n / 2, (n / 5) / (2 * _c2), (n / 5) / (2 * 2)
    v0 = voigt(0, g, l)
    yv = [a * voigt(n - p0, g, l) / v0  for n in range(0, n)]
    coeff = [1, - 1 / n,  0 / n ** 2]
    yp = [poly(el, *coeff) for el in range(n)]
    yn = np.random.randn(n) * std
    y = np.sum([yv, yp, yn], axis = 0)
    yf = fit_voigt(range(n), y, peak_ratio = 1/5, background = 1)[-1]


# Add Docs
fit_voigt.__doc__  = fit_voigt_doc
fit_gaussian.__doc__  = fit_gaussian_doc
fit_lorentian.__doc__ = fit_lorentian_doc