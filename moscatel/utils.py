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

        if skip_every is not None:
            print('Skipping every {0}-th frames per band\n'.format(skip_every))

        else: #elif skip_every == None:
            '''
            does not print even if skip_every is not entered in terminal
            '''
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
        '''
        bug: title shows as band_idx instead of g,r,z
        '''
        axx = df[cols].plot(subplots=True, alpha=0.8, figsize=(15,8))
        plt.suptitle('Raw flux in {}-band'.format(band_idx))
        #axx.set_ylabel('Raw Flux')
        #axx.set_xlabel('Time (HJD)')

    plt.show()
    return df

def plot_multicolor(df, star_idx, showfig=None):

    if star_idx==0:
        cols = 'g_a_flux g_b_flux g_c_flux'.split()

    elif star_idx==1:
        cols = 'r_a_flux r_b_flux r_c_flux'.split()

    else:
        cols = 'z_a_flux z_b_flux z_c_flux'.split()

    if showfig is not None and showfig == True:#showfig==None or showfig==True:
        #normalize
        df = df/df.median()
        import pdb; pdb.set_trace()
        #df.index.to_julian_date().values
        axs = df[cols].plot(subplots=False, color=['g','r','b'], figsize=(15,8), marker='o')
        axs.set_ylabel('Raw Flux')
        axs.set_xlabel('Time (HJD)')

    plt.show()


def df_phot(target, ref, df, normed, showfig=None):
    '''
    small bug: imperfect normalization
    either caused by dataframe conversion or median does not correspond to baseline
    '''
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
    #convert df to series to fix imperfect normalization problem
    res=(df[t]/df[r]).values

    #normalization
    if normed == True:
        res /= np.median(res)
    #res = pd.DataFrame({t.column: res}, index=df.index)
    else:
        pass
    if showfig==None or showfig==True:
        fig, ax2 = plt.subplots(nrows=1,ncols=1,figsize=(10,5))
        plt.plot(df.index, res, 'ko', linestyle='none')
        ax2.set_title('{0}-band of {1}/{2}'.format(df.columns[0].split('_')[0],
        target,ref))
        # ax2 = res.plot(figsize=(15,5), color='k', marker='o', linestyle='none', title='{}-band'.format(df.columns[0].split('_')[0]))
        ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%m'))
        ax2.set_ylabel('Normalized Flux')
        ax2.set_xlabel('Time (HJD)')
        plt.show()
    return df[t]/df[r]

def plot_matrix(df):
    scatter_matrix(df, figsize=(10,10), marker='o', alpha=0.5);
    plt.show()

def df_phot_multicolor(target, ref, df_g, df_r, df_z, normed, showfig=None):
    res1 = df_phot(target=target, ref=ref, df=df_g, normed=normed, showfig=False)
    res2 = df_phot(target=target, ref=ref, df=df_r, normed=normed, showfig=False)
    res3 = df_phot(target=target, ref=ref, df=df_z, normed=normed, showfig=False)
    #convert series to dataframe
    res1 =res1.to_frame()
    res2 =res2.to_frame()
    res3 =res3.to_frame()

    df_grz = res1.join([res2, res3])
    if showfig == True:
        #ax3 = df_grz.iloc[:-20].plot(figsize=(15,5), marker='o', legend=False, linestyle='none', title='g-,r-,z-band of {0}/{1}'.format(target,ref))
        ax3 = df_grz.plot(figsize=(10,5), marker='o', color=['g','r','b'], alpha=0.8, legend=False, linestyle='none', title='g-,r-,z-band of {0}/{1}'.format(target,ref))
        patches, labels = ax3.get_legend_handles_labels()
        ax3.legend(patches, ['g','r','z'], loc='best')
        ax3.xaxis.set_major_formatter(dates.DateFormatter('%H:%m'))
        ax3.set_ylabel('Normalized Flux')
        ax3.set_xlabel('Time (HJD)')
        plt.show()
    return df_grz

def df_phot_multicolor2(target, ref, df_g, df_r, df_z, star, normed, showfig=None):
    res1 = df_phot(target=target, ref=ref, df=df_g, normed=normed, showfig=True)
    res2 = df_phot(target=target, ref=ref, df=df_r, normed=normed, showfig=True)
    res3 = df_phot(target=target, ref=ref, df=df_z, normed=normed, showfig=True)
    return res1, res2, res3

# def df_phot_multicolor2(target, ref, df_grz, star, normed, showfig=None):
#     if target=='a':
#         t=df_grz.columns[[0,10,18]]
#     elif target=='b':
#         t=df_grz.columns[[3,12,21]]
#     else:
#         t=df_grz.columns[[6,15,24]]
#     if ref=='a':
#         r=df_grz.columns[[0,10,18]]
#     elif ref=='b':
#         r=df_grz.columns[[3,12,21]]
#     else:
#         r=df_grz.columns[[6,15,24]]
#     #differential photometry
#     res_grz=df_grz[t]/df_grz[r]
#
#     #normalization
#     if normed == True:
#         #(res-res.mean())/res.std()
#         res_grz = res_grz/res_grz.median().astype(np.float64)
#     else:
#         pass
#     if showfig==None or showfig==True:
#         plot_multicolor(res_grz, star, showfig=True)
#     return res_grz

def plot_params(df):
    print('Parameters to choose from:\n{}'.format(df.columns))
    try:
        cols = raw_input('columns to plot (e.g. z_a_x z_a_y): ')
    except:
        cols = input('columns to plot (e.g. z_a_x z_a_y): ')
    cols = cols.split()
    df[cols].plot(subplots=True, figsize=(15,4))
    plt.show()

def check_data(output_dir):
    print('\nNOTE:\n')
    try:
        df_g = pd.read_csv(output_dir+'/gband.csv', index_col=0, parse_dates=True)
    except:
        print(output_dir+'/gband.csv does not exist')
    try:
        df_r = pd.read_csv(output_dir+'/rband.csv', index_col=0, parse_dates=True)
    except:
        print(output_dir+'/rband.csv does not exist')
    try:
        df_z = pd.read_csv(output_dir+'/zband.csv', index_col=0, parse_dates=True)
    except:
        print(output_dir+'/zband.csv does not exist')
    # try:
    #     df_grz = pd.read_csv(output_dir+'/grzband.csv', index_col=0, parse_dates=True)
    # except:
    #     print(output_dir+'/grzband.csv does not exist')

    return df_g, df_r, df_z


def remove_outlier(df_g, df_r, df_z):
    print('\n-----------------------')
    print('Removing outliers outside +/-sigma={} in each band...'.format(clip_sigma))
    print('-----------------------')
    df_g = df_g[np.abs(df_g-df_g.mean())<=(clip_sigma*df_g.std())]
    n_g = df_g[~(np.abs(df_g-df_g.mean())<=(clip_sigma*df_g.std()))]
    print('removed {} outliers in g-band'.format(len(n_g)))
    df_r = df_r[np.abs(df_r-df_r.mean())<=(clip_sigma*df_r.std())]
    n_r = df_r[~(np.abs(df_r-df_r.mean())<=(clip_sigma*df_r.std()))]
    print('removed {} outliers in r-band'.format(len(n_r)))
    n_z = df_z[~(np.abs(df_z-df_z.mean())<=(clip_sigma*df_z.std()))]
    df_z = df_z[np.abs(df_z-df_z.mean())<=(clip_sigma*df_z.std())]
    print('removed {} outliers in z-band'.format(len(n_z)))
    #check if grz
    #df_grz = df_grz[np.abs(df_grz-df_grz.mean())<=(clip_sigma*df_grz.std())]
    return df_g, df_r, df_z
