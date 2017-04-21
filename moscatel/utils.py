#!/usr/bin/env python

import numpy as np
from photutils.centroids import centroid_com as com
from photutils import CircularAperture
from photutils import CircularAnnulus
from photutils import aperture_photometry
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as dates
from pandas.tools.plotting import scatter_matrix

import os
from glob import glob
from tqdm import tqdm

try:
	from astropy.io import fits as pf
except:
	import pyfits as pf

def init_moscatel(filedir, skip_every=None):
	file_list = glob(os.path.join(filedir,'*.fits'))
	file_list.sort()
	if os.listdir(filedir) != []:
	#if len(file_list)>0:
		print('total no. of raw data frames: {0}\n'.format(len(file_list)))

		if skip_every:
			print('Skipping {0}th frames raw frames per band'.format(skip_every))

		elif skip_every == None:
			print('Analyzing all raw frames per band')
			skip_every=1

		gband=[]
		rband=[]
		zband=[]

		#get list of frames by filter based on header
		for i in tqdm(file_list[::skip_every]):
			hdr = pf.open(i)[0].header
			if hdr['FILTER'] == 'g':
				gband.append(i)
			elif hdr['FILTER'] == 'r':
				rband.append(i)
			else: #hdr['FILTER'] == 'z_s':
				zband.append(i)

		print('gband={0} frames\nrband={1} frames\nzband={2} frames'.format(len(gband), len(rband), len(zband)))
		bands=(gband,rband,zband)

	else:
		print('ERROR: check your data directory')
		bands=''

	return bands

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

def plot_lightcurve(dfs, band_idx, showfig=None):
    df = pd.concat(dfs, axis=1)
    #df.head()

    if band_idx==0:
        cols = 'g_a_flux g_b_flux g_c_flux'.split()

    elif band_idx==1:
        cols = 'r_a_flux r_b_flux r_c_flux'.split()

    else:
        cols = 'z_a_flux z_b_flux z_c_flux'.split()

    if showfig==None or showfig==True:
        axx = df[cols].plot(subplots=True, figsize=(15,8))
        #axx.set_ylabel('Raw Flux')
        #axx.set_xlabel('Time (HJD)')

    plt.show()
    return df

def plot_multicolor(df, star_idx, showfig=None):

    if star_idx==0:
        cols = 'g_a_flux r_a_flux z_a_flux'.split()

    elif star_idx==1:
        cols = 'g_b_flux r_b_flux z_b_flux'.split()

    else:
        cols = 'z_c_flux z_c_flux z_c_flux'.split()

    if showfig==None or showfig==True:
    	#normalize
    	df = df/df.max().astype(np.float64)
    	axs = df[cols].plot(subplots=False,figsize=(15,8), marker='o')
        #axs.set_ylabel('Raw Flux')
        #axs.set_xlabel('Time (HJD)')

    plt.show()


def df_phot(target, ref, df, normed, showfig=None):
    if target=='a':
        t=df.columns[0]
    elif target=='b':
        t=df.columns[3]
    else:
        t=df.columns[6]

    if ref=='a':
        r=df.columns[0]
    elif ref=='b':
        r=df.columns[3]
    else:
        r=df.columns[6]

    #differential photometry
    res=df[t]/df[r]

    #normalization
    if normed == True:
        #(res-res.mean())/res.std()
        res = res/res.max().astype(np.float64)
    else:
        pass
    if showfig==None or showfig==True:
        ax2 = res.plot(figsize=(15,5), color='k', marker='o', linestyle='none', title='{}-band'.format(df.columns[0].split('_')[0]))
        ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%m'))
        ax2.set_ylabel('Normalized Flux')
        ax2.set_xlabel('Time (HJD)')
        plt.show()
    return res

def plot_matrix(df):
    scatter_matrix(df, figsize=(15,15), marker='o', alpha=0.5);
    plt.show()

def df_phot_multicolor(target, ref, df_grz, star, normed, showfig=None):
    if target=='a':
        t=df_grz.columns[[0,10,18]]
    elif target=='b':
        t=df_grz.columns[[3,12,21]]
    else:
        t=df_grz.columns[[6,15,24]]

    if ref=='a':
        r=df_grz.columns[[0,10,18]]
    elif ref=='b':
        r=df_grz.columns[[3,12,21]]
    else:
        r=df_grz.columns[[6,15,24]]

    #differential photometry
    res_grz=df_grz[t]/df_grz[r]

    #normalization
    if normed == True:
        #(res-res.mean())/res.std()
        res_grz = res_grz/res_grz.max().astype(np.float64)
    else:
        pass
    if showfig==None or showfig==True:
        plot_multicolor(res_grz, star, showfig=True)

    return res_grz

def plot_params(df):
    print('Parameters to choose from:\n{}'.format(df.columns))
    try:
        cols = raw_input('columns to plot (e.g. z_a_x z_a_y): ')
    except:
        cols = input('columns to plot (e.g. z_a_x z_a_y): ')
    cols = cols.split()
    df[cols].plot(subplots=True, figsize=(15,4))
    plt.show()
