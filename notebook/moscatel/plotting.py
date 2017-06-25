<<<<<<< HEAD
import matplotlib.pyplot as plt
import pandas as pd
from photutils import CircularAperture
from astropy.visualization import ZScaleInterval
=======
import pandas as pd
from photutils import CircularAperture
from astropy.visualization import ZScaleInterval
import warnings
from moscatel import utils
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
import matplotlib.pyplot as plt
>>>>>>> d1e37f14af9fbee428799658187fe419867b8809
import numpy as np

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
<<<<<<< HEAD
=======

def show_one_image_all_centroid(img, centroid, key, nstars=10, figsize=(10,5)):
    #create new figure every band
    fig = plt.figure(figsize=figsize)
    #centroid = list(zip(sources[band_id]['xcentroid'],sources[band_id]['ycentroid']))
    for i,xy in enumerate(centroid[:nstars]):
        #bkg subtraction?
                
        try:
            img_crop = utils.get_crop(img, xy, box_size=40)
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

def show_one_centroid_all_image(band, centroid, skip_every, ncols=8, figsize=(10,8)):
    
    nrows=round(len(band[::skip_every])/ncols)
    fig = plt.figure(figsize=figsize)
    for idx,i in enumerate(band[::skip_every]):
        img=pf.open(i)[0].data
        hdr=pf.open(i)[0].header
        try:
            img_crop = utils.get_crop(img, centroid, box_size=60)
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
>>>>>>> d1e37f14af9fbee428799658187fe419867b8809
