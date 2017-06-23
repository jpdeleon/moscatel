#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from tqdm import tqdm
import pandas as pd
import getpass
import ast
from astropy.visualization import ZScaleInterval

try:
    from astropy.io import fits as pf
except:
    import pyfits as pf

def init_moscatel(filedir, filters_in_config, output_dir, skip_every=None):
    """Checks the raw data and sorts them by band/color.

        Parameters
        ----------

        filedir : str
                path to raw data directory

        skip_every : int
                interval of raw data to skip for quick look;
                e.g., skip_every=2 skips every 2nd raw data;
                default is 1
    """
    file_list = glob(os.path.join(filedir,'*.fits'))
    file_list.sort()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if os.listdir(filedir) != []:
    #if len(file_list)>0:
        print('total no. of raw data frames: {0}\n'.format(len(file_list)))

        if skip_every is not None:
            print('Skipping every {0}-th frames per band\n'.format(skip_every))

        else: #elif skip_every == None:
            '''
            bug: does not print even if skip_every is not entered in terminal
            '''
            print('Analyzing all raw frames per band')
            skip_every=1

        bands = {}
        for j in filters_in_config:
            j=j.strip(' ')
            bands[j]=[]
        filters_in_hdr=[]

        #get list of frames by filter based on header
        for i in tqdm(file_list[::skip_every]):
            hdr = pf.open(i)[0].header
            filters_in_hdr.append(hdr['FILTER'])
            for j in filters_in_config:
                if hdr['FILTER'] == j:
                    j=j.strip(' ')
                    bands[j].append(i)

        filters_in_hdr_set=list(set(filters_in_hdr)).sort()

        for k in bands.keys():
            print('{0}-band={1} frames'.format(k, len(bands[k])))
            #save into txtfile
            name = os.path.join(output_dir,k+'-band.txt')
            with open(name, 'w') as z: #overwrite
                #write line by line
                for line in bands[k]:
                    z.write('{}\n'.format(line))
        print('\nfilenames sorted by band saved in {}'.format(output_dir))

    else:
        print('ERROR: check your data directory')
        #return empty dict
        bands={}

    return filters_in_hdr_set, bands

def create_config(config_dir):
    '''
    Called in moscatel-init to create
    a default configuration file (if does not exist)
    '''
    fname=os.path.join(config_dir,'config.txt')
    with open(fname, 'w') as w: #overwrite
        w.write('#moscatel configuration file\n')
        w.write('#folders must be at home/user/\n')
        w.write('#remove any whitespace\n')
        w.write('#-------------------------------------------------------\n')
        w.write('#initialize: moscatel-init\n')
        w.write('data_dir=data/hatp44_data\n')
        w.write('output_dir=output\n')
        w.write('filters=g,r,z_s\n')
        w.write('#-------------------------------------------------------\n')
        w.write('#photometry settings: moscatel-phot\n')
        w.write('centroids=(703, 303),(915, 264),(707, 758)\n')
        w.write('aperture_radius=20\n')
        w.write('#-------------------------------------------------------\n')
        w.write('#lightcurve analysis settings: moscatel-analysis\n')
    print('config.txt created in {}'.format(config_dir))

def check_config():
    home_dir=os.path.join('/home',getpass.getuser(),'moscatel')
    config=np.genfromtxt(os.path.join(home_dir,'config.txt'),dtype=str, \
                        delimiter='=', comments='#')
    for i in config:
        if i[0] == 'data_dir':
            data_dir = os.path.join('/home',getpass.getuser(),i[-1].strip(' '))
        elif i[0] == 'output_dir':
            output_dir = os.path.join(data_dir,i[-1].strip(' '))
        elif i[0] == 'filters':
            filters = i[-1].strip(' ')
        elif i[0] == 'centroids':
            #convert txt to tuple
            centroids = ast.literal_eval(i[-1])
        elif i[0] == 'aperture_radius':
            aper_radius = i[-1]
            aper_radius = ast.literal_eval(aper_radius)
        else:
            print('config file not read correctly')
            data_dir, output_dir, filters, centroids, aper_radius = [],[],[],[],[]
    return data_dir, output_dir, filters, centroids, aper_radius


def combine_df(output_dir, clip):
    df_r = pd.read_csv(output_dir+'/rband.csv', index_col=0, parse_dates=True)
    df_g = pd.read_csv(output_dir+'/gband.csv', index_col=0, parse_dates=True)
    df_z = pd.read_csv(output_dir+'/zband.csv', index_col=0, parse_dates=True)
    #combine df and save
    df_grz = df_g.join([df_r, df_z])
    filename=output_dir+'/grzband.csv'
    df_grz.to_csv(filename, mode = 'w', header =df_grz.columns)
    #read
    df_grz = pd.read_csv(output_dir+'/grzband.csv', index_col=0, parse_dates=True)
    #df.index.to_julian_date()
    if clip is not None:
        try:
            '''
            uneven clipping happens because len of each dataset is different;
            e.g., g-band becomes shorter than r-band lc
            '''
            df_g = df_g.iloc[clip[0]:-clip[1]]
            df_r = df_r.iloc[clip[0]:-clip[1]]
            df_z = df_z.iloc[clip[0]:-clip[1]]
            df_grz = df_grz.iloc[clip[0]:-clip[1]]
        except:
            print('ERROR: cannot clip: {0},{1}'.format(clip[0], clip[1]))
    return df_r, df_g, df_z, df_grz

def load_df(output_dir, band, clip):
    try:
        fname='{}-band_phot.csv'.format(band)
        df = pd.read_csv(os.path.join(output_dir,fname), index_col=0, parse_dates=True)
        if clip is not None:
            df = df.iloc[clip[0]:-clip[1]]
    except:
        print('\nNOTE: check missing {0}-band data in {1}'.format(band, output_dir))

    return df

def clip_outlier(clip_sigma, t, r, df_g, df_r, df_z, df_grz):
    '''
    bug: number of outliers (in target,ref cols) not correctly counted
    e.g. try: n_g = df_g['g_'+t+'_flux'][...]
    '''
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
    try:
        df_grz = df_grz[np.abs(df_grz-df_grz.mean())<=(clip_sigma*df_grz.std())]
    except:
        pass

    return df_g, df_r, df_z, df_grz
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

def stack_raw_image(image_list, skip_every=1):
    '''
    stack image using median to be used for detecting
    source locations (i.e. target and ref stars)
    '''
    image_array = []
    for i in tqdm(image_list[::skip_every]):
        img = pf.getdata(i)
        image_array.append(img)
    stacked_image = np.median(image_array, axis=0)

    # stacked_image_g = stack_raw_image(bands['g'], skip_every=10)
    # stacked_image_r = stack_raw_image(bands['r'], skip_every=10)
    # stacked_image_z = stack_raw_image(bands['z_s'], skip_every=10)

    return stacked_image

def show_stacked_images(images):
    '''

    '''
    fig, axes = plt.subplots(1,len(images),figsize=(15,5))
    titles=list(images.keys())
    for i,band in enumerate(images):
        vmin,vmax= ZScaleInterval().get_limits(images[band])
        axes[i].imshow(images[band],vmin=vmin,vmax=vmax)
        axes[i].set_title(titles[i])
    plt.show()
    #return None

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

def save_df(input_dir,df,band_names,band_idx):
    #save dataframe as csv
    filename=os.path.join(input_dir,
    '{0}-band_phot.csv'.format(band_names[band_idx]))

    if os.path.isfile(filename):
        print('\nOverwriting {}'.format(filename))


    df.to_csv(filename, mode = 'w', header =df.columns)

    print('\n-----------------------')
    print('Data saved in {}'.format(filename))
    print('-----------------------\n')
