#!/usr/bin/env python

from utils import get_crop, get_centroid, get_phot, get_bkg, get_phot2
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

def make_lightcurve(centroids, bands, band_idx, box_size):
    for star_idx in range(3):
        xcenters, ycenters = [],[]
        aperture_sums = []
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

            #compute background
            #bkg_mean=get_bkg(image_crop, centroid, r_in=20., r_out=30.)

            #without aperture photometry
            aperture_sum = get_phot(image_crop, centroid, r=20)

            #minus background wihtin annulus
            #aperture_sum = get_phot2(image_crop,bkg_mean,centroid,r=20)

            aperture_sums.append(aperture_sum)

        #output as dataframe of given band and star
        dfs.append(pd.DataFrame(
            {'{0}_{1}_x'.format(band_name[band_idx], star_names[star_idx]) : xcenters,
             '{0}_{1}_y'.format(band_name[band_idx], star_names[star_idx]) : ycenters,
             '{0}_{1}_flux'.format(band_name[band_idx], star_names[star_idx]) : aperture_sums},
            index = obs_time))

    return dfs, band_idx
