import warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit
from scipy.ndimage.filters import gaussian_filter
from moscatel import utils

def gauss1D(x, *params):
    A, mu, sigma, eps= params
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + eps

def model_gaussian(image_crop,convolve=False,verbose=False, show_fit=False, return_fwhm=False):
    #normalize image
    if convolve==True:
        sigma_estimate = 4
        image_crop = gaussian_filter(image_crop,sigma=sigma_estimate)
    image_crop /= np.max(image_crop)
    # https://python4astronomers.github.io/fitting/sherpa.html
    i,j = np.unravel_index(image_crop.argmax(), image_crop.shape) #take x,y max
    peak_x=image_crop[i,:]
    peak_y=image_crop[:,j]
    
    #estimate mean and standard deviation
    ydata = (peak_x+peak_y)/2.0
    #ydata /= np.max(ydata)
    xdata = np.array(range(len(ydata)))
    xmean = len(xdata)/2.0
    sigma = np.std(ydata)
    amp = np.max(ydata)
    eps = np.median(ydata)
    #import pdb;pdb.set_trace()
    #fitting
    popt, pcov = curve_fit(gauss1D, xdata, ydata, p0 = [amp, xmean, sigma, eps])

    if show_fit == True:
        plt.plot(xdata,gauss1D(xdata, *popt), label='Gaussian fit')
        plt.plot(xdata,ydata,'o',label='data',alpha=0.5)
        plt.legend()
    if verbose==True:
        print('A: {}\nmu: {}\nsigma= {}\neps: {}'.format(popt[0],popt[1], popt[2], popt[3]))

    if return_fwhm==True:
        fwhm=utils.sigma_to_fwhm(popt[2])
        return popt, fwhm
    return popt


def model_gaussian2D(img_crop, verbose=False, fwhm=8., return_fwhm=False):
    sigma= utils.fwhm_to_sigma(fwhm)
    try:
        #get 1D fit results
        result_1Dfit = model_gaussian(img_crop)
        amp, mu, sigma, eps = result_1Dfit
        #initialize model
        g_init = models.Gaussian2D(amplitude=amp,x_mean=mu, y_mean=mu, x_stddev=sigma, y_stddev=sigma)
    except:
	sigma= utils.fwhm_to_sigma(fwhm)
        #if 1D fitting fails due to wrong initial centroiding
        x_mean, y_mean = img_crop.shape[0]/2, img_crop.shape[1]/2 
        #initialize model using default values
        g_init = models.Gaussian2D(amplitude=1,x_mean=x_mean, y_mean=y_mean, x_stddev=sigma, y_stddev=sigma)

    #fit model
    fit_g = fitting.LevMarLSQFitter()

    #normalize image before fitting
    img_crop_norm =img_crop/np.max(img_crop)

    x,y=range(img_crop.shape[0]),range(img_crop.shape[1])
    with warnings.catch_warnings():
        # Ignore model linearity warning from the fitter
        warnings.simplefilter('ignore')
        g = fit_g(g_init, x, y, img_crop_norm[y,x])

    if verbose==True:
        print(g.param_names,g.param_sets)
        
    if return_fwhm==True:
        fwhm_mean = utils.sigma_to_fwhm(np.abs((g.x_stddev.value+g.y_stddev.value)/2))
        return g, fwhm_mean	

    return g

def get_fwhm(img_crop, convolve=False, verbose=False, show_fit=False, method='1D'):
    try:
        fit_result1D = model_gaussian(img_crop,convolve=convolve, verbose=verbose, show_fit=show_fit);
    except:
        warnings.warn('no good 1D gaussian fit')
    if method=='2D':
        try:
            fit_result2D = model_gaussian2D(img_crop,convolve=convolve, verbose=verbose, show_fit=show_fit);
        except:
            warnings.warn('no good 2D gaussian fit')
        g = model_gaussian2D(img_crop)
        sigma_mean = np.abs((g.x_stddev.value+g.y_stddev.value)/2)
        fwhm_mean = utils.sigma_to_fwhm(sigma_mean)
        return fwhm_mean
    
    else: #1D
        return utils.sigma_to_fwhm(fit_result1D[2])
