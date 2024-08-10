The package **`peaky`** allows the user to fit a single peak to a Gaussian, Lorentian or Voigt profile. 

## Basic Example
Here is a basic way to fit your data to a Voigt profile:
```
import peaky
peaky.fit_voigt(x, y)
```
here `x` and `y` are respectively the independent and dependent variable. The data should contain a single peak, not larger then the data set. 
If instead you would like to specifically fit the data to a gaussian or a lorentian profile, use one of these two functions instead:
```
peaky.fit_gaussian(x, y)
peaky.fit_lorentian(x, y)
```

## Parameters
These are all other parameters of the `fit_voigt` function:

- **background**
It sets the order of the polynomial to helps the fit adapt to a possible background signal. A value of 0 means constant background, a value of 1 means linear background etc. If None is provided, no background is added. The dafault value is 1. 

- **a**
It sets the amplitude of the voigt peak to the value provided. In this context, amplitude means literally the amplitude of the signal at peak position. If None is provided (as by default) the best amplitude is found by the fit. 

- **p**
It sets the peak position to the value provided. If None is provided (as by default) the best peak position is found by the fit. 

- **g**
It sets the voigt gaussian width to the value provided. If None is provided (as by default) the best gaussian width is found by the fit. 

- **l**
It sets the voigt lorentian width to the value provided. If None is provided (as by default) the best lorentian width is found by the fit. 

- **peak_ratio**
This parameters helps the process find the best initial values for the fit. For example, a value of 0.3 means then the peak occupies roughly 1 / 3 of the data set width (estimated by eye). If a lower value is provided, the noise may ruin the peak detection, while a higher value may average away the signal.

- **log**
When `True` (as by default), all measured parameters are printed directly on terminal together with their error bars and guessed initial values. A color coded fit signal to noise ratio may help decide if a fit is reliable. Red means unreliable, orange means not optimal, while green means reliable.

## Output
These are the outputs of the fit functions:

- **fit**
The optimal values of the fit parameters in this order: amplitude, peak position, gaussian and lorentian width and the full width at half maximum (fwhm). To these, it follows the list of polynomial coefficients used for the background signal from lower to higher order. Note that each background coefficient (being very small for higher orders) is multiplied by the n ** i where n is the length of the signal and i the order of the coefficient.

- **err**
The error correspondent to each parameters obtained from the fit.

- **snr**
The fit signal to noise measured as the ratio between the signal amplitude (obtained from fit) and the fit data noise (signal - fit standard deviation). This parameter may help decide whatever or not a fit is reliable.

- **fit functions**
A list of three functions: the first corresponding to the overall fitted signal, the second only to the Voigt component (without background) and the third to the background signal alone.

- **fit data**
The correspondent data values (three lists) correspondent to the three previous functions evaluated in x, which could be easily used in a plot.

## Gaussian and Lorentian Fit
Note that the `fit_gaussian` function is equivalent to `fit_voigt` with the Lorentzian width parameter `l` fixed to 0, while `fit_lorentian` is equivalent to the `fit_voigt` with Gaussian width parameter `g` fixed to 0. 

## Example
Here is the plot of a simulated data set (a = 1, p = 1500, g = 150, l = 100 with linear background) and of the correspondent fit performed by peaky:

![example](https://github.com/piccolomo/peaky/raw/master/images/plot.png)

And here is the printed output:

![example](https://github.com/piccolomo/peaky/raw/master/images/output.png)

## Test
You can run a simple test of the newly installed package, to check that `peaky` works well in your machine. Just use `peaky.run_test()`

## Other Documentation
The documentation of the functions shown above could be accessed using the following commands:
```
print(peaky.fit_voigt.__doc__)
print(peaky.fit_gaussian.__doc__)
print(peaky.fit_lorentian.__doc__)
```

## Installation
To install the latest version of the `peaky` package use the following command:
```
sudo -H pip3 install peaky
```

## Credits
- Author: Savino Piccolomo
- Inpired during a research visit at the Laboratory of Nano Optics in Siegen
- e-mail: piccolomo@gmail.com
