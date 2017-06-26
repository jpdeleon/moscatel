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
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
from moscatel import utils
from moscatel import models

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
    centroid_com(): Calculates the object “center of mass” from 2D image moments.
    centroid_1dg(): Calculates the centroid by fitting 1D Gaussians to the marginal x and y distributions of the data.
    centroid_2dg(): Calculates the centroid by fitting a 2D Gaussian to the 2D distribution of the data.
    Default is centroid_2dg.
    centroid_com(): Calculates the object center of mass from 2D image moments
    centroid_1dg(): Calculates the centroid by fitting 1D Gaussians to the marginal x and y distributions of the data
    centroid_2dg(): Calculates the centroid by fitting a 2D Gaussian to the 2D distribution of the data
    Default is centroid_2dg
    ''' 
    if method=='com':
        x, y = centroid_com(image)
    
    elif method=='1d_gaussian':
        x, y = centroid_1dg(image)

    else: #default
        x, y = centroid_2dg(image)
        
    return (x,y)

def get_bkg(image, centroid, r_in, r_out):
    '''
    To calculate the mean local background within the circular annulus aperture
    with thickness dr, divide its sum by its area calculated using the area() method
    '''
    annulus = CircularAnnulus(centroid, r_in, r_out)
    result = aperture_photometry(image, annulus)
    bkg_mean = float(result['aperture_sum'] / annulus.area())
    '''
    The background sum within the circular aperture is then the mean local background 
    times the circular aperture area
    
    bkg_sum = bkg_mean * apertures.area()
    residual=aperture_sum - bkg_sum
    '''
    return bkg_mean #, residual


def get_phot(image, centroid, r=10):
    '''
    r is based from measured fwhm = 8.0
    '''
    #if isinstance(r, int):
    #elif isinstance(r, tuple):
    apertures = CircularAperture(centroid, r)
    phot_table = aperture_photometry(image, apertures)
    
    #xcenter = phot_table['xcenter']
    #ycenter = phot_table['ycenter']
    #centroid = (xcenter, ycenter)
    aperture_sum = float(phot_table['aperture_sum'])
    
    return aperture_sum #,centroid

def get_phot_table(image, centroid, aperture_radii):#, error):
    '''
    aperture_sum_err provides the propagated uncertainty associated with 'aperture_sum
    '''    
    apertures = [CircularAperture(centroid, r=r) for r in aperture_radii]
    
    err=0.01*image
    
    phot_table = aperture_photometry(image, apertures, method='subpixel', error=err)
    #convert to pandas dataframe
    '''
    bug: Cannot convert a table with mixin columns to a pandas DataFrame
    '''
    #phot_table = phot_table.to_pandas()
    return phot_table

def get_phot2(image, bkg_mean, centroid, r=10):
        
    apertures = CircularAperture(centroid, r)
    phot_table = aperture_photometry(image - bkg_mean, apertures)
    aperture_sum = float(phot_table['aperture_sum'])
    
    return aperture_sum #,centroid


def get_annulus(img_crop, r_in=15, r_out=20, show_mask=False):
    xcenter,ycenter = img_crop.shape[0]/2, img_crop.shape[1]/2

    y, x = np.indices((img_crop.shape))
    mask1 = np.sqrt((x - xcenter)**2 + (y - ycenter)**2) <= r_out
    mask2 = np.sqrt((x - xcenter)**2 + (y - ycenter)**2) <= r_in
    if show_mask:
        plt.imshow(mask1^mask2) #subtract boolean
    return (mask1^mask2)

def get_uncertainty(img_crop, r_in=15, r_out=20, mode='std'):
    mask=get_annulus(img_crop, r_in, r_out)
    if mode=='std':
        return np.std(img_crop[mask])
    else:
        print('unknown mode')
        return None

def get_peak_flux(img_crop, r_in=15, r_out=20):
    mask=get_annulus(img_crop, r_in, r_out)
    return np.max(img_crop[mask])

def make_lightcurve(band_list, star_positions, aperture_radii, skip_every=1, box_size=40, r_in=15, r_out=20):#, save_as_df=False):
    xcenters, ycenters = [], []
    uncertainty, peak_flux = [], []
    tables = {}
    frames={}
    print('performing aperture photometry on {} stars\n'.format(len(star_positions)))
    for star_idx, position in enumerate(star_positions):
        tables[star_idx]=[]
        print('\n---------star index: {}---------'.format(star_idx))
        print('initial centroid: {}'.format(position))
        #each star position determined from the stacked image
        for i in tqdm(band_list[::skip_every]):
            #each image in a given band
            image = pf.getdata(i)
            hdr = pf.getheader(i)
            #crop
            img_crop = utils.get_crop(image, position, box_size=box_size)
            #compute new centroid
            new_centroid = get_centroid(img_crop, method='2D_gaussian')
            #perform photometry
            phot_table = get_phot_table(img_crop, new_centroid, aperture_radii=aperture_radii)
            #estimate background within annulus
            bkg= get_bkg(img_crop, new_centroid, r_in, r_out)
            phot_table['bkg'] = bkg
            #add uncertainty within prefefined annulus
            unc= get_uncertainty(img_crop, r_in=15, r_out=20, mode='std')
            phot_table['uncertainty'] = unc
            #get peak flux within the same annulus
            peak_flux=get_peak_flux(img_crop, r_in=15, r_out=20)
            phot_table['peak_flux'] = peak_flux
            '''
            BUG: image_crop is normalized due to rescaling in get_fwhm() computation 
            '''
            #measure fwhm
            fwhm = models.get_fwhm(img_crop)
            phot_table['fwhm'] = fwhm
            #calculate centroid shift 
            xshift= phot_table['xcenter'].value - box_size/2.
            yshift= phot_table['ycenter'].value - box_size/2.
            #correct actual centroid (not cropped image)
            phot_table['xcenter'] = xshift+position[0]
            phot_table['ycenter'] = yshift+position[1]
            
            xcenters.append(phot_table['xcenter'])
            ycenters.append(phot_table['ycenter'])
            #add time (as index)
            phot_table['mjd'] = hdr['MJD-STRT']
            tables[star_idx].append(phot_table)
            #import pdb; pdb.set_trace()

    print('---------Done---------')
    return tables


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
