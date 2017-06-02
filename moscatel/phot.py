#!/usr/bin/env python
from photutils.centroids import centroid_com as com
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
import numpy as np
from tqdm import tqdm
try:
    from astropy.io import fits
except:
    import pyfits as fits
from datetime import datetime as dt
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

def get_crop(image, centroid, box_size):
    x, y = centroid
    image_crop = np.copy(image[int(y-(box_size/2)):int(y+(box_size/2)),int(x-(box_size/2)):int(x+(box_size/2))])

    return image_crop

def get_centroid(image):
    '''
    Calculate the centroid of a 2D array as its 'center of mass' determined from image moments.
    '''
    centroid = com(image)
    return centroid

def get_phot(image, centroid, r):
    fwhm = 8.0

    apertures = CircularAperture(centroid, r)
    phot_table = aperture_photometry(image, apertures)

    #xcenter = phot_table['xcenter']
    #ycenter = phot_table['ycenter']
    #centroid = (xcenter, ycenter)
    aperture_sum = float(phot_table['aperture_sum'])

    return aperture_sum #,centroid

def sigma_per_r(image, centroid, r_in, r_out, delta_r,show_image=True):
    r = np.arange(r_in,r_out,delta_r)
    aperture_sums = []
    for i in r:
        aperture_sums.append(get_phot(image, centroid, r=i))
    if show_image==True:
        plt.plot(r,aperture_sums,'o')
        plt.xlabel('aperture radius')
        plt.ylabel('aperture sum')
    return aperture_sums

def radial_profile(image, center):
    y, x = np.indices((image.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), image.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile

def get_bkg(image, centroid, r_in=10., r_out=20.):
    annulus = CircularAnnulus(centroid, r_in, r_out)
    result = aperture_photometry(image, annulus)
    bkg_mean = float(result['aperture_sum'] / annulus.area())
    return bkg_mean

def get_phot2(image, bkg_mean, centroid, r=10):

    apertures = CircularAperture(centroid, r)
    phot_table = aperture_photometry(image - bkg_mean, apertures)
    aperture_sum = float(phot_table['aperture_sum'])

    return aperture_sum #,centroid

dfs = []
def make_lightcurve(centroids, bands, band_idx, box_size, aperture_radius):
    """Creates a lightcurve after doing aperture photometry

        Parameters
        ----------

        centroids : typle
                (x,y) centroid positions from config.txt

        bands: dict
                dictionary of filenames by band/color

        band_idx: int
                index for band/color

        box_size: int
                size of box for cropping (in pix)

        aperture_radius: int
                aperture radius for photometry (in pix)
    """
    band_names = np.sort(list(bands.keys()))
    num_stars= range(len(centroids))
    for star_idx in num_stars:
        xcenters, ycenters = [],[]
        aperture_sums = []
        background = []
        fwhms = []
        obs_time = []
        obs_mjd = []
        ##extract lightcurve (enumerate all frames) in a given band
        for i in tqdm(bands[band_names[band_idx]]):
            #import pdb; pdb.set_trace()
            hdr = fits.open(i)[0].header
            img = fits.open(i)[0].data
            #get dates from fits header
            date=dt.strptime(hdr['DATE-OBS'], '%Y-%m-%d')
            time=dt.strptime(hdr['EXP-STRT'], '%H:%M:%S.%f')
            newdate = time.replace(year=date.year, month=date.month, day=date.day)
            obs_time.append(newdate)
            obs_mjd.append(hdr['MJD-STRT'])

            #crop
            #import pdb; pdb.set_trace()
            image_crop = get_crop(img, centroids[star_idx], box_size)

            ###aperture photometry###
            #compute centroid
            centroid = get_centroid(image_crop)

            xcenters.append(centroid[0])
            ycenters.append(centroid[1])

            #compute backgound
            bkg_mean=get_bkg(image_crop, centroid, r_in=20., r_out=30.)

            #measure fwhm
            fwhm=get_fwhm(image_crop)

            #without aperture photometry

            aperture_sum = get_phot(image_crop, centroid, r=aperture_radius)

            #minus background wihtin annulus
            #aperture_sum = get_phot2(image_crop,bkg_mean,centroid,r=aperture_radius)

            aperture_sums.append(aperture_sum)
            background.append(bkg_mean)

            # if fwhm < 10*np.median(fwhms):
            #     fwhms.append(fwhm)
            # else:
            #     fwhms.append(np.nan)
            fwhms.append(fwhm)

        #output as dataframe of given band and star

        dfs.append(pd.DataFrame(
            {'{0}_{1}_x'.format(band_names[band_idx], str(star_idx)) : xcenters,
             '{0}_{1}_y'.format(band_names[band_idx], str(star_idx)) : ycenters,
             '{0}_{1}_flux_r{2}'.format(band_names[band_idx], str(star_idx), aperture_radius) : aperture_sums,
             '{0}_{1}_bkg_r{2}'.format(band_names[band_idx], str(star_idx), aperture_radius) : background,
             '{0}_{1}_fwhm_r{2}'.format(band_names[band_idx], str(star_idx), aperture_radius) : fwhms},
             #'airmass' : airmass
            index = obs_time))
    return dfs, band_idx, band_names

def get_fwhm(image_crop):
    # https://python4astronomers.github.io/fitting/sherpa.html
    i,j = np.unravel_index(image_crop.argmax(), image_crop.shape) #take x,y max
    #plt.plot(image_crop[i,:], '-')
    #plt.plot(image_crop[:,j], '--')
    peak_x=image_crop[i,:]
    peak_y=image_crop[:,j]
    try:
        sigma=model_gaussian(peak_x, peak_y)
        fwhm=2.355*np.abs(sigma)
    except:
        #no good estimate
        fwhm=np.nan

    return fwhm

def gauss(x, *params):
    A, mu, sigma, eps= params
    return A*np.exp(-(x-mu)**2/(2.*sigma**2)) + eps

def model_gaussian(peak_x, peak_y,verbose=False):
    #estimate mean and standard deviation
    ydata = (peak_x+peak_y)/2.0
    xdata = np.array(range(len(ydata)))
    mean = np.mean(ydata)
    sigma = np.std(ydata)
    amp = np.max(ydata)
    eps =0.1
    #fitting
    popt, pcov = curve_fit(gauss, xdata, ydata, p0 = [amp, mean, sigma, eps])

    #plt.plot(xdata,gauss(xdata, *popt), label='Gaussian fit')
    #plt.plot(xdata,ydata,'ok', label='data')
    #plt.legend()
    if verbose==True:
        print('A: {}\nmu: {}\nsigma= {}\neps: {}'.format(popt[0],popt[1], popt[2], popt[3]))
    return popt[2]
