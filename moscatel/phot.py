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

dfs = []
band_name = ['g','r','z']
star_names = 'abc'

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
    bkg_mean = result['aperture_sum'] / annulus.area()
    return bkg_mean

def get_phot2(image, bkg_mean, centroid, r=10):

    apertures = CircularAperture(centroid, r)
    phot_table = aperture_photometry(image - bkg_mean, apertures)
    aperture_sum = float(phot_table['aperture_sum'])

    return aperture_sum #,centroid

def make_lightcurve(centroids, bands, band_idx, box_size, aperture_radius):
    for star_idx in range(3):
        xcenters, ycenters = [],[]
        aperture_sums = []
        background = []
        obs_time = []
        obs_mjd = []
        sum_per_band = {}

        ##extract lightcurve (enumerate all frames) in a given band
        for i in tqdm(bands[band_idx]):
            hdr = fits.open(i)[0].header
            img = fits.open(i)[0].data

            #get dates from fits header
            date=dt.strptime(hdr['DATE-OBS'], '%Y-%m-%d')
            time=dt.strptime(hdr['EXP-STRT'], '%H:%M:%S.%f')
            newdate = time.replace(year=date.year, month=date.month, day=date.day)
            obs_time.append(newdate)
            obs_mjd.append(hdr['MJD-STRT'])

            #crop
            image_crop = get_crop(img, centroids[star_idx], box_size)

            ###aperture photometry###
            #compute centroid
            centroid = get_centroid(image_crop)
            centroids.append(centroid)

            xcenters.append(centroid[0])
            ycenters.append(centroid[1])

            #compute backgound
            bkg_mean=get_bkg(image_crop, centroid, r_in=20., r_out=30.)

            #without aperture photometry
            aperture_sum = get_phot(image_crop, centroid, r=aperture_radius)

            #minus background wihtin annulus
            #aperture_sum = get_phot2(image_crop,bkg_mean,centroid,r=aperture_radius)

            aperture_sums.append(aperture_sum)
            background.append(bkg_mean)

        #output as dataframe of given band and star
        '''
        fix column to add aperture radius: e.g.
        '{0}_{1}_flux_r{2}'.format(band_name[band_idx], star_names[star_idx], aperture_radius) : aperture_sums},
        '''
        dfs.append(pd.DataFrame(
            {'{0}_{1}_x'.format(band_name[band_idx], star_names[star_idx]) : xcenters,
             '{0}_{1}_y'.format(band_name[band_idx], star_names[star_idx]) : ycenters,
             '{0}_{1}_flux'.format(band_name[band_idx], star_names[star_idx]) : aperture_sums},
            index = obs_time))

    return dfs, band_idx
