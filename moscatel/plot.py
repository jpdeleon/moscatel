#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime as dt
import matplotlib.dates as dates
import numpy as np
from pandas.tools.plotting import scatter_matrix

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


def df_phot(target, ref, df, normed=True, showfig=None):
    '''
    small bug: imperfect normalization
    either caused by dataframe conversion or median does not correspond to baseline

    revise structure of df: arbitrary number of ref stars
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

def plot_params(df):
    print('Parameters to choose from:\n{}'.format(df.columns))
    try:
        cols = raw_input('columns to plot (e.g. z_a_x z_a_y): ')
    except:
        cols = input('columns to plot (e.g. z_a_x z_a_y): ')
    cols = cols.split()
    df[cols].plot(subplots=True, figsize=(15,4))
    plt.show()

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
