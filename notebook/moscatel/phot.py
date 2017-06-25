from photutils import DAOStarFinder
from photutils.detection import find_peaks
from photutils import CircularAperture
from astropy.stats import sigma_clipped_stats
from astropy.visualization import ZScaleInterval
from photutils import centroid_com, centroid_1dg, centroid_2dg
from photutils import CircularAnnulus
from photutils import aperture_photometry
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

def get_sources(img, num_stars=10, fwhm=8.0):
    mean, median, std = sigma_clipped_stats(img, sigma=3.0, iters=5)
    daofind = DAOStarFinder(fwhm=fwhm, threshold=5.*std)  
    sources = daofind(img - median)
    #convert to pandas dataframe for easy sorting
    sources = sources.to_pandas()
    #sort by brightness
    sources = sources.sort_values(by='peak',ascending=False)
    
    return sources.head(num_stars)


def get_centroid(image, method='com'):
    '''
<<<<<<< HEAD
    centroid_com(): Calculates the object center of mass from 2D image moments
    centroid_1dg(): Calculates the centroid by fitting 1D Gaussians to the marginal x and y distributions of the data
    centroid_2dg(): Calculates the centroid by fitting a 2D Gaussian to the 2D distribution of the data
    Default is centroid_2dg
=======
    centroid_com(): Calculates the object “center of mass” from 2D image moments.
    centroid_1dg(): Calculates the centroid by fitting 1D Gaussians to the marginal x and y distributions of the data.
    centroid_2dg(): Calculates the centroid by fitting a 2D Gaussian to the 2D distribution of the data.
    Default is centroid_2dg.
>>>>>>> 9135c664d4486dfdf5e12b8717f6ea2b27239bd4
    ''' 
    if method=='com':
        x, y = centroid_com(image)
    
    elif method=='1d_gaussian':
        x, y = centroid_1dg(image)

    else: #default
        x, y = centroid_2dg(image)
        
    return (x,y)

def get_bkg(image, centroid, r_in=10., r_out=20.):
    annulus = CircularAnnulus(centroid, r_in, r_out)
    result = aperture_photometry(image_crop, annulus)
    bkg_mean = result['aperture_sum'] / annulus.area()
    return bkg_mean


def get_phot(image, centroid, r=10):
    '''
    r is based from measured fwhm = 8.0
    '''
    
    apertures = CircularAperture(centroid, r)
    phot_table = aperture_photometry(image, apertures)
    
    #xcenter = phot_table['xcenter']
    #ycenter = phot_table['ycenter']
    #centroid = (xcenter, ycenter)
    aperture_sum = float(phot_table['aperture_sum'])
    
    return aperture_sum #,centroid

def get_phot2(image, bkg_mean, centroid, r=10):
        
    apertures = CircularAperture(centroid, r)
    phot_table = aperture_photometry(image - bkg_mean, apertures)
    aperture_sum = float(phot_table['aperture_sum'])
    
    return aperture_sum #,centroid


def df_phot(target, ref, df, showfig=None):
    if target=='a':
        t=df.columns[0]
    elif target=='b':
        t=df.columns[3]
    else:
        t=df.columns[6]
    if ref=='a':
        r=df.columns[0]
    elif target=='b':
        r=df.columns[3]
    else:
        r=df.columns[6]
        
    res=df[t]/df[r]
    
    fig, ax2 = plt.subplots(1,1,figsize=(10,8))
    if showfig==None or showfig==True:
        res.plot(figsize=(15,5), color='k', marker='o', linestyle='none', ax=ax2);
        
    return res
