import numpy as np
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
import matplotlib.pyplot as plt

import pandas as pd

from photutils import CircularAperture
from astropy.visualization import ZScaleInterval
import warnings
from moscatel import phot
from moscatel import utils
from moscatel import models
from astropy.stats import sigma_clipped_stats
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage.filters import gaussian_filter

def show_sources(image, sources, labels, method='sources', num_stars=10):
    '''
    similar to `show_peaks`; difference is `sources` used as input
    '''
    if isinstance(image, np.ndarray):
        if method == 'peaks':
            positions = (peaks['x_peak'].values[:num_stars], peaks['y_peak'].values[:num_stars])
        else: #default
            positions = (sources['xcentroid'].values[:num_stars], sources['ycentroid'].values[:num_stars])
        
        apertures = CircularAperture(positions, r=20.)
        vmin,vmax= ZScaleInterval().get_limits(image)
        plt.figure(figsize=(10,10))
        plt.imshow(image, origin='lower', vmin=vmin,vmax=vmax)
        for num, (x,y) in enumerate(zip(positions[0],positions[1])):
            plt.text(x+5,y+5, num+1, fontsize=20, color='w')
        apertures.plot(color='r', lw=2)
            
    elif len(image) == len(sources) and isinstance(image, list):
        fig, ax = plt.subplots(1,len(image),figsize=(15,5))
        for idx,(img,src) in enumerate(zip(image,sources)):
            positions = (src['xcentroid'].values[:num_stars], src['ycentroid'].values[:num_stars])
            apertures = CircularAperture(positions, r=20.)
            vmin,vmax= ZScaleInterval().get_limits(img)
            ax[idx].imshow(img, origin='lower', vmin=vmin,vmax=vmax)
            for num, (x,y) in enumerate(zip(positions[0],positions[1])):
                ax[idx].text(x+5,y+5, num+1, fontsize=20, color='w')
            apertures.plot(color='r', lw=2, ax=ax[idx])
        for idx, label in enumerate(labels):
            ax[idx].set_title(label)
            
    else:
        print('incorrect dimensions')
    plt.show()
    #return None

def plot_lightcurve(dfs, band_idx, showfig=None):
    df = pd.concat(dfs, axis=1)
    #df.head()

    if band_idx==0:
        cols = 'g_a_flux g_b_flux g_c_flux'.split()

    elif band_idx==1:
        cols = 'r_a_flux r_b_flux r_c_flux'.split()

    else:
        cols = 'z_a_flux z_b_flux z_c_flux'.split()
    
    fig, ax = plt.subplots(1,1,figsize=(15,10))
    if showfig==None or showfig==True:
        df[cols].plot(subplots=True, figsize=(15,8),ax=ax)
        
    return df

def show_one_image_all_centroid(img, centroid, key, nstars=10, boxsize=40, figsize=(10,5)):
    #create new figure every band
    fig = plt.figure(figsize=figsize)
    #centroid = list(zip(sources[band_id]['xcentroid'],sources[band_id]['ycentroid']))
    for i,xy in enumerate(centroid[:nstars]):
        #bkg subtraction?
                
        try:
            img_crop = utils.get_crop(img, xy, box_size=boxsize)
            ncols=nstars/2
            nrows=int(nstars/ncols+1)
            ax = plt.subplot(nrows,ncols,i+1)
        except:
            warnings.warn('star not cropped')
        ax.set_title('star idx={}'.format(i))
        ax.imshow(img_crop)
        ax.plot(img_crop.shape[0]/2,img_crop.shape[1]/2, 'r+', markersize=10)
        ax.set_xlabel('{0:.1f},  {1:.1f}'.format(xy[0],xy[1]))
        ax.set_xticks([])
        ax.set_yticks([])
    plt.suptitle('{}-band'.format(key))
    #return None

def show_one_centroid_all_image(band, centroid, skip_every, boxsize=40, ncols=8, figsize=(10,8)):
    
    nrows=round(len(band[::skip_every])/ncols)
    fig = plt.figure(figsize=figsize)
    for idx,i in enumerate(band[::skip_every]):
        img=pf.open(i)[0].data
        hdr=pf.open(i)[0].header
        try:
            img_crop = utils.get_crop(img, centroid, box_size=boxsize)
            ax = plt.subplot(nrows,ncols,idx+1)
        except:
            warnings.warn('star not cropped')
        ax.set_title('img={}'.format(idx))
        ax.imshow(img_crop)
        ax.plot(img_crop.shape[0]/2,img_crop.shape[1]/2, 'w+', markersize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    plt.suptitle('{}-band'.format(hdr['FILTER']))
    plt.axis('tight')
    plt.show()
    #return None 


def show_fit_2D(sample_img, centroid, sigma_estimate=4, boxsize=40, show_image=True, convolve=False, recenter=False, show_3D=False):
    if convolve==True:
        sample_img = gaussian_filter(sample_img, sigma=sigma_estimate)

    if recenter == True:
        #crop a bigger box first
        img_crop = utils.get_crop(sample_img, centroid, box_size=80)
        #############RECENTROID#############
        xy_new = phot.get_centroid(img_crop, method='2D_gaussian')

        #------------RE-CROP with smaller box------------#
        img_crop = utils.get_crop(img_crop, xy_new, box_size=boxsize)
    else:
        img_crop = utils.get_crop(sample_img,centroid,box_size=boxsize)

    fig = plt.figure(figsize=(15,5))

    #define and fit model
    g = models.model_gaussian2D(img_crop)
    xx,yy=np.mgrid[0:img_crop.shape[0],0:img_crop.shape[1]]

    if show_image==True:
        #data
        ax1 = plt.subplot(1,3,1)
        ax1.imshow(img_crop, origin='lower', interpolation='nearest')
        ax1.contour(img_crop, colors='w')
        ax1.plot(img_crop.shape[0]/2,img_crop.shape[1]/2, 'r+', markersize=20)
        ax1.set_title('data')
        ax1.set_xlabel('{0:.1f},  {1:.1f}'.format(centroid[0],centroid[1])) #centroid x,y
        #model
        ax2 = plt.subplot(1,3,2)
        ax2.imshow(g(yy,xx), origin='lower', interpolation='nearest')
        ax2.contour(img_crop, colors='w')
        ax2.plot(img_crop.shape[0]/2,img_crop.shape[1]/2, 'r+', markersize=20)
        ax2.set_title('model')
        ax2.set_xlabel(g.param_sets[1:3].flatten()) #centroid x,y
        #residual
        ax3 = plt.subplot(1,3,3)
        ax3.imshow(img_crop-g(yy,xx), origin='lower', interpolation='nearest')
        ax3.plot(img_crop.shape[0]/2,img_crop.shape[1]/2, 'r+', markersize=20)
        ax3.set_title('residual')

        '''
        Take note of the 3D syntax: g(yy,xx)[yy,xx]
        '''
        if show_3D==True:
            fig = plt.figure(figsize=(15,5))
            #data
            img_crop_norm = img_crop/np.max(img_crop)
            ax1 = fig.add_subplot(131, projection='3d')
            ax1.plot_surface(X=xx, Y=yy, Z=img_crop_norm[yy,xx])
            ax1.set_title('data (rescaled)')
            #model
            ax2 = fig.add_subplot(132, projection='3d')
            ax2.plot_surface(X=xx, Y=yy, Z=g(yy,xx)[yy,xx])
            ax2.set_title('model')
            #residual
            ax3 = fig.add_subplot(133, projection='3d')
            ax3.plot_surface(X=xx, Y=yy, Z=(img_crop-g(yy,xx)[yy,xx]))
            ax3.set_title('residual')
    return g

def plot_psf(img_crop):
    mean, median, stddev=sigma_clipped_stats(img_crop)
    vmin,vmax= ZScaleInterval().get_limits(img_crop)

    fig = plt.figure(figsize=(15,5))
    ax0 = fig.add_subplot(131)
    im=ax0.imshow(img_crop-median,vmin=vmin,vmax=vmax)
    fig.colorbar(im,ax=ax0)
    
    xdim=img_crop.shape[0]/2
    xslice=img_crop[int(xdim),:]
    ax1 = fig.add_subplot(132)
    ax1.plot(xslice/np.max(xslice))

    ax2 = fig.add_subplot(133, projection='3d')
    #ax = plt.gca(projection='3d')
    xx,yy=np.mgrid[0:img_crop.shape[0],0:img_crop.shape[1]]

    ax2.plot_surface(X=xx, Y=yy, Z=img_crop[yy,xx])
    #return None
