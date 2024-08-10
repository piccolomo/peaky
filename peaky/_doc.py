fit_voigt_doc = """It fits the data to a single voigt profile with polynomial background. The data should contain a single peak, not larger then the data set. Here is a basic example:

   \x1b[92mimport peaky\x1b[0m
   \x1b[92mpeaky.fit_voigt(x, y)\x1b[0m

where x and y are respectively the independent and dependent variable. Here are all other parameters: 

\x1b[94mbackground\x1b[0m
It sets the order of the polynomial to helps the fit adapt to a possible background signal. A value of 0 means constant background, a value of 1 means linear background etc. If None is provided, no background is added. The dafault value is 1. 

\x1b[94ma\x1b[0m
It sets the amplitude of the voigt peak to the value provided. In this context, amplitude means literally the amplitude of the signal at peak position. If None is provided (as by default) the best amplitude is found by the fit. 

\x1b[94mp\x1b[0m
It sets the peak position to the value provided. If None is provided (as by default) the best peak position is found by the fit. 

\x1b[94mg\x1b[0m
It sets the voigt gaussian width to the value provided. If None is provided (as by default) the best gaussian width is found by the fit. 

\x1b[94ml\x1b[0m
It sets the voigt lorentian width to the value provided. If None is provided (as by default) the best lorentian width is found by the fit. 

\x1b[94mpeak_ratio\x1b[0m
This parameters helps the process find the best initial values for the fit. For example, a value of 0.3 means then the peak occupies roughly 1 / 3 of the data set width (estimated by eye). If a lower value is provided, the noise may ruin the peak detection, while a higher value may average away the signal. 

\x1b[94mlog\x1b[0m
When True (as by default), all measured parameters are printed directly on terminal together with their error bars and guessed initial values. A color coded fit signal to noise ratio may help decide if a fit is reliable. Red means unreliable, orange means not optimal, while green means reliable. 

\nThese are the outputs of the fit:

\x1b[33mfit\x1b[0m
The optimal values of the fit parameters in this order: amplitude, peak position, gaussian and lorentian width and the full width at half maximum (fwhm). To these, it follows the list of polynomial coefficients used for the background signal from lower to higher order. Note that each background coefficient (being very small for higher orders) is multiplied by the n ** i where, n is the length of the signal and i is the order of the coefficient. 

\x1b[33merr\x1b[0m
The errors for each of the previous parameters obtained from the fit.

\x1b[33msnr\x1b[0m
The fit signal to noise measured as the ratio between the signal amplitude (obtained from fit) and the fit data noise. This parameter may help decide whatever or not a fit is reliable. 

\x1b[33myf\x1b[0m
The list of fitted data values that can be used, for example, for a plot.

\nAll measured parameters are printed directly on terminal together with their error bars and guessed initial values. A color coded fit signal to noise ratio may help decide if a fit is reliable. Red means unreliable, orange means not optimal, while green means reliable. 

"""

fit_gaussian_doc = """It is equivalent to the fit_voigt function with lorentian width parameter l fixed to 0. Access the fit_voigt function doc string for full documentation. """

fit_lorentian_doc = """It is equivalent to the fit_voigt function with gaussian width parameter g fixed to 0. Access the fit_voigt function doc string for full documentation. """
