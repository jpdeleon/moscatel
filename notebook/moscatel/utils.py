import numpy as np
from tqdm import tqdm
try:
    from astropy.io import fits as pf
except:
    import pyfits as pf
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval

#suppose these come form config file
filters_in_config = 'g,r,z_s'.split(',')
bands = {}
def get_band_list(file_list):
    for j in filters_in_config:
        #initialize dict with empty arrays
        j=j.strip(' ')
        bands[j]=[]
    filters_in_hdr=[]

    for i in tqdm(file_list):
        hdr = pf.getheader(i)
        filters_in_hdr.append(hdr['FILTER'])
        for j in filters_in_config:
            if hdr['FILTER'] == j:
                j=j.strip(' ')
                bands[j].append(i)

    for key in sorted(bands.keys()):
        print('{0}-band: {1} frames'.format(key, len(bands[key])))
    return bands


def stack_raw_images(image_list, skip_every=1):
    '''
    stack image using median to be used for detecting 
    source locations (i.e. target and ref stars)
    '''
    image_array = []
    count = []
    for i in tqdm(image_list[::skip_every]):
        img = pf.getdata(i)
        image_array.append(img)
        count.append(i)
    stacked_image = np.median(image_array, axis=0)
    print('number of stacked raw images={}'.format(len(count)))
    return stacked_image

def show_stacked_images(images):
    fig, axes = plt.subplots(1,3,figsize=(15,5))
    titles='g,r,z'.split(',')
    for i,img in enumerate(images):
        vmin,vmax= ZScaleInterval().get_limits(img)
        axes[i].imshow(img,vmin=vmin,vmax=vmax)
        axes[i].set_title(titles[i])
    plt.show()
    #return None

def get_crop(image, centroid, box_size):
    x, y = centroid
    image_crop = np.copy(image[int(y-(box_size/2)):int(y+(box_size/2)),int(x-(box_size/2)):int(x+(box_size/2))])

    return image_crop


def fwhm_to_sigma(fwhm):
    return fwhm/ (2*np.sqrt(2*np.log(2)))

def sigma_to_fwhm(sigma):
    return sigma * (2*np.sqrt(2*np.log(2)))

if __name__ == "__main__":
    main()
