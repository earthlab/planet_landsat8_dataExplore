import os, sys
import rasterio
from rasterio.mask import mask
import numpy as np
from matplotlib import pyplot as plt
import glob
from xml.dom import minidom
import fiona
from shapely.geometry import mapping, shape
import scipy.misc as misc
import matplotlib.lines as mlines
from mpl_toolkits.axes_grid1 import make_axes_locatable

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.15)
    return fig.colorbar(mappable, cax=cax)

# function to return collection dates from Planet MultiSpectral image metadata
def getDates_planet(xml_file):

    xmldoc = minidom.parse(xml_file)
    date_node = xmldoc.getElementsByTagName("eop:acquisitionDate")
    date_str = str(date_node[0].firstChild.nodeValue)
    acq_date_yyyymmdd = date_str.split('T')[0]        
    
    return acq_date_yyyymmdd

# function to return correction coefficients from Planet MultiSpectral image metadata
def getCorrCoefs_planet(xml_file):

    xmldoc = minidom.parse(xml_file)
    nodes = xmldoc.getElementsByTagName("ps:bandSpecificMetadata")

    # XML parser refers to bands by numbers 1-4
    coeffs = {}
    for node in nodes:
        bn = node.getElementsByTagName("ps:bandNumber")[0].firstChild.data
        if bn in ['1', '2', '3', '4']:
            i = int(bn)
            value = node.getElementsByTagName("ps:reflectanceCoefficient")[0].firstChild.data
            coeffs[i] = float(value)
            
    return coeffs

# function to extract image patch
def getImagePatch_planet(imfile, poly):
    
    try:
        with rasterio.open(imfile, 'r') as src:
            arr,_ = mask(src, [poly], crop=True)

        if arr.sum() != 0.:
            return arr
        else:
            return 0
        
    except:
        return 0

# function return geometries as geoJSON
def getGeometries_planet(shpfile):
    
    with fiona.open(shpfile, "r") as shapefile:
        geoms = [feature["geometry"] for feature in shapefile]
        
    return geoms


def getMeansBand(arr, nbands):
    
    res=[]
    for i in range(nbands):
        res.append(arr[i,:,:].mean())
    
    return res

def getSTDband(arr, nbands):
    
    res=[]
    for i in range(nbands):
        res.append(arr[i,:,:].std())
    
    return res

def cor_ms(ms_temp, cf_arr):
    # extract the band means per patch
    ms_band1 = [b[0,:,:] for b in ms_temp]
    ms_band2 = [b[1,:,:] for b in ms_temp]
    ms_band3 = [b[2,:,:] for b in ms_temp]
    ms_band4 = [b[3,:,:] for b in ms_temp]
    
    # correct the ms images to TOA reflectance
    ms_band1_cor = [cf[0] * mn for cf,mn in zip(cf_arr, ms_band1)]
    ms_band2_cor = [cf[1] * mn for cf,mn in zip(cf_arr, ms_band2)]
    ms_band3_cor = [cf[2] * mn for cf,mn in zip(cf_arr, ms_band3)]
    ms_band4_cor = [cf[3] * mn for cf,mn in zip(cf_arr, ms_band4)]
    
    res = []
    for i in range(len(ms_band1_cor)):
        res.append(np.rollaxis(np.dstack((ms_band1_cor[i], ms_band2_cor[i], ms_band3_cor[i], ms_band4_cor[i])),2,0))
    
    return res
    
def meanIntensityOverGeometry_simple_Landsat(bandlist, poly):
    
    # SR scale factor
    sr_scale = 10000.
    
    # get the band patches. should be 4 of size 1xrowsxcols items, [blue, green, red, nir]
    band_patches = [getImagePatch_planet(im, poly) for im in bandlist]
    #band_temp = [np.ma.masked_equal(arr, -9999.)/sr_scale for arr in band_patches if type(arr) is np.ma.core.MaskedArray]
    band_temp = [np.ma.masked_equal(arr, -9999.)/sr_scale for arr in band_patches]
    #band_means = [getMeansBand(im,1) for im in band_temp]
    band_means = [im.mean() for im in band_temp]
    band_std = [im.std() for im in band_temp]
    
    res = {}
    res['sr_means'] = band_means
    res['numpixels'] = band_temp[0].count()
    res['sr_std'] = band_std
    
    #band_temp = [np.ma.masked_equal(arr, -9999.) for arr in band_temp]
    res['image'] = np.ma.array([b.squeeze() for b in band_temp])
    return res
    
# define function which process the planet images within a specified geometry
def meanIntensityOverGeometry(home_dir, poly):
    
    all_ims = glob.glob('{}/*/*.tif'.format(home_dir))
#     all_ms_ims = [im for im in all_ims if "MS.tif" in im] # only want ms tifs associated with SR files
    all_sr_ims = [im for im in all_ims if "SR" in im]
    all_qa = [im for im in all_ims if "DN_UDM" in im]
#     all_xml = glob.glob('{}/*/*.xml'.format(home_dir)) # only want the xmls associated with SR files
    all_xml = [s.replace('MS_SR.tif', 'MS_metadata.xml') for s in all_sr_ims]
    all_ms_ims = [s.replace('MS_SR.tif', 'MS.tif') for s in all_sr_ims]
    
    #print(all_xml)
    
    c = getCorrCoefs_planet(all_xml[0])
    
    ms_patches = [getImagePatch_planet(im, poly) for im in all_ms_ims]
    sr_patches = [getImagePatch_planet(im, poly) for im in all_sr_ims]
    qa_patches = [getImagePatch_planet(im, poly) for im in all_qa]
    coeffs = [getCorrCoefs_planet(xml) for xml in all_xml]
    dates = [getDates_planet(xml) for xml in all_xml]
    
    #print(dates)
    
    # extract the mean values per band. however, we may be getting nodata (0) pixels influencing the mean.
    # may need to break this out into separate steps
    xml_inds = [i for i,j in enumerate(sr_patches) if type(j) is np.ma.core.MaskedArray]
    xml_inds = [i for i,j in enumerate(sr_patches) if np.sum(j) > 10]
    
    #print(xml_inds)
    
    # put the correction coefficients into a usable list corresponding to only overlapping images
    cf_arr=[]
    for i in xml_inds:
        cf_arr.append(np.array([coeffs[i][1], coeffs[i][2], coeffs[i][3], coeffs[i][4]]))
    
    

    ms_patch_temp = [ms_patches[i] for i in xml_inds]
    ms_temp = [np.ma.masked_equal(arr, 0.0) for arr in ms_patch_temp if type(arr) is np.ma.core.MaskedArray]
    ms_cor = cor_ms(ms_temp, cf_arr) # correct the MS imagery to TOA reflectance
    ms_means = [getMeansBand(arr,4) for arr in ms_cor]
    ms_means = [getMeansBand(arr,4) for arr in ms_patches]

    sr_patch_temp = [sr_patches[i] for i in xml_inds]
    sr_temp = [np.ma.masked_equal(arr, 0.0) for arr in sr_patch_temp if type(arr) is np.ma.core.MaskedArray]
    sr_means = [getMeansBand(arr,4) for arr in sr_temp]
    sr_means = [getMeansBand(arr,4) for arr in sr_patches]
    
    # extract the band means per patch
    ms_band1_means = [b[0] for b in ms_means]
    ms_band2_means = [b[1] for b in ms_means]
    ms_band3_means = [b[2] for b in ms_means]
    ms_band4_means = [b[3] for b in ms_means]

#     print sr_means[0]
    sr_scale = 10000
    sr_band1_means = [b[0]/sr_scale for b in sr_means]
    sr_band2_means = [b[1]/sr_scale for b in sr_means]
    sr_band3_means = [b[2]/sr_scale for b in sr_means]
    sr_band4_means = [b[3]/sr_scale for b in sr_means]
    
    # put the dates in a usable list corresponding to only overlapping images
    date_arr = [dates[i] for i in xml_inds]
    
    #print(1)

    res = {}
    res['ms_means'] = [ms_band1_means, ms_band2_means, ms_band3_means, ms_band4_means]
    res['sr_means'] = [sr_band1_means, sr_band2_means, sr_band3_means, sr_band4_means]
    res['dates'] = date_arr
    res['sr_image'] = sr_patches
    res['ms_image'] = ms_patches
    res['qa_image'] = qa_patches
    
#     res['ms_std'] = [ms_band1_std_cor, ms_band2_std_cor, ms_band3_std_cor, ms_band4_std_cor]
#     res['sr_std'] = [sr_band1_std, sr_band2_std, sr_band3_std, sr_band4_std]
    return res
    
# define function which process the planet images within a specified geometry
def meanIntensityOverGeometry_simple(im, xml, poly):
    
    # SR scale factor
    sr_scale = 10000.    
    
    #ms_patch = getImagePatch_planet(im, poly)
    sr_patch = getImagePatch_planet(im, poly) #/ sr_scale 
    
    coeffs = getCorrCoefs_planet(xml)
    date_s = getDates_planet(xml)
    
    #ms_patch_temp = [ms_patches[i] for i in xml_inds]
    #ms_temp = [np.ma.masked_equal(arr, 0.0) for arr in ms_patch_temp if type(arr) is np.ma.core.MaskedArray]
    #ms_cor = cor_ms(ms_temp, cf_arr) # correct the MS imagery to TOA reflectance
    #ms_means = [getMeansBand(arr,4) for arr in ms_cor]

    #sr_patch_temp = [sr_patches[i] for i in xml_inds]
    #sr_temp = [np.ma.masked_equal(arr, 0.0) for arr in sr_patch_temp if type(arr) is np.ma.core.MaskedArray]
    sr_temp = np.ma.masked_equal(sr_patch, 0.0) / sr_scale
    sr_means = getMeansBand(sr_temp,4)
    sr_stds = getSTDband(sr_temp, 4)
    
    
    
    # extract the band means per patch
    #ms_band1_mean = [b[0] for b in ms_means]
    #ms_band2_mean = [b[1] for b in ms_means]
    #ms_band3_mean = [b[2] for b in ms_means]
    #ms_band4_mean = [b[3] for b in ms_means]

    
    sr_band1_mean = sr_means[0] #/sr_scale 
    sr_band2_mean = sr_means[1] #/sr_scale 
    sr_band3_mean = sr_means[2] #/sr_scale 
    sr_band4_mean = sr_means[3] #/sr_scale
    
    sr_band1_std = sr_stds[0]
    sr_band2_std = sr_stds[1]
    sr_band3_std = sr_stds[2]
    sr_band4_std = sr_stds[3]
    
    # put the dates in a usable list corresponding to only overlapping images
    #date_arr = [dates[i] for i in xml_inds]

    res = {}
    #res['ms_means'] = [ms_band1_means, ms_band2_means, ms_band3_means, ms_band4_means]
    res['sr_means'] = [sr_band1_mean, sr_band2_mean, sr_band3_mean, sr_band4_mean]
    res['sr_std'] = [sr_band1_std, sr_band2_std, sr_band3_std, sr_band4_std]
    res['dates'] = date_s
    res['numpixels'] = sr_temp[0,:,:].count()
    res['image'] = sr_temp
#     res['ms_std'] = [ms_band1_std_cor, ms_band2_std_cor, ms_band3_std_cor, ms_band4_std_cor]
#     res['sr_std'] = [sr_band1_std, sr_band2_std, sr_band3_std, sr_band4_std]
    return res


def returnSamples(pl_imArr, oli_imArr, method='direct'):
    '''this function returns the planet samples within a landsat pixel.
        The list of methods is to account for how the landsat pixel is sampled for each set of planet pixels.
        Currently, 'direct' is the only one and assigns the same landsat pixel for each planet pixel in the set'''
    
    # get the rows and columns of the landsat data
    oli_rows, oli_cols = oli_imArr.shape
    
    # get the rows and columns of the planet data
    planet_rows, planet_cols = pl_imArr.shape
    
    # get the row/column factors for how many planet pixels will fit in a landsat pixel
    row_fact = int(np.floor(planet_rows/oli_rows))
    col_fact = int(np.floor(planet_cols/oli_cols))
    
    # it is not necessarily an even number given the spatial extraction, so add in the modulus
    row_mod = planet_rows%oli_rows
    col_mod = planet_cols%oli_cols
    
    # get the ranges of pixels to grab from planet
    row_sets = np.arange(0,planet_rows, row_fact) + row_mod
    col_sets = np.arange(0, planet_cols, col_fact) + col_mod
    
    # sample via landsat
    samples = []
    for r in range(oli_rows-1):
        for c in range(oli_cols-1):
            
            # get the landsat pixel
            oli_pix_val = oli_imArr[r,c]
            
            # get the pixels from planet if the pixel is not 0 or -9999
            if oli_pix_val is not np.ma.masked: # set the condition later
                pl_row_start = row_sets[r]
                pl_row_end = row_sets[r+1]
                pl_col_start = col_sets[c]
                pl_col_end = col_sets[c+1]
                pl_pix = pl_imArr[pl_row_start:pl_row_end, pl_col_start:pl_col_end].ravel()
                
                oli_pix = np.ones(pl_pix.shape)*oli_pix_val
                samples.append(np.vstack((oli_pix, pl_pix)))
    
    return samples

def returnSamplesMeanStd(pl_imArr, oli_imArr, method='direct'):
    '''this function returns the planet samples within a landsat pixel.
        The list of methods is to account for how the landsat pixel is sampled for each set of planet pixels.
        Currently, 'direct' is the only one and assigns the same landsat pixel for each planet pixel in the set'''
    
    # get the rows and columns of the landsat data
    oli_rows, oli_cols = oli_imArr.shape
    
    # get the rows and columns of the planet data
    planet_rows, planet_cols = pl_imArr.shape
    
    # get the row/column factors for how many planet pixels will fit in a landsat pixel
    row_fact = int(np.floor(planet_rows/oli_rows))
    col_fact = int(np.floor(planet_cols/oli_cols))
    
    # it is not necessarily an even number given the spatial extraction, so add in the modulus
    row_mod = planet_rows%oli_rows
    col_mod = planet_cols%oli_cols
    
    # get the ranges of pixels to grab from planet
    row_sets = np.arange(0,planet_rows, row_fact) + row_mod
    col_sets = np.arange(0, planet_cols, col_fact) + col_mod
    
    # sample via landsat
    samples = []
    for r in range(oli_rows-1):
        for c in range(oli_cols-1):
            
            # get the landsat pixel
            oli_pix_val = oli_imArr[r,c]
            
            # get the pixels from planet if the pixel is not 0 or -9999
            if oli_pix_val is not np.ma.masked: # set the condition later
                pl_row_start = row_sets[r]
                pl_row_end = row_sets[r+1]
                pl_col_start = col_sets[c]
                pl_col_end = col_sets[c+1]
                pl_pix = pl_imArr[pl_row_start:pl_row_end, pl_col_start:pl_col_end].ravel()
                pl_pix_mean = np.mean(pl_pix)
                pl_pix_std = np.std(pl_pix)
                pl_pix_max = np.max(pl_pix)
                pl_pix_min = np.min(pl_pix)
                
                
                samples.append(np.array([[oli_pix_val], [pl_pix_mean], [pl_pix_std], [pl_pix_max - pl_pix_min]]))
    
    return samples

def returnSamplesMeanStd_ndvi(pl_imArr_red, pl_imArr_nir, oli_imArr_red, oli_imArr_nir, method='direct'):
    '''this function returns the planet samples within a landsat pixel.
        The list of methods is to account for how the landsat pixel is sampled for each set of planet pixels.
        Currently, 'direct' is the only one and assigns the same landsat pixel for each planet pixel in the set'''
    
    # get the rows and columns of the landsat data
    oli_rows, oli_cols = oli_imArr_red.shape
    
    # get the rows and columns of the planet data
    planet_rows, planet_cols = pl_imArr_red.shape
    
    # get the row/column factors for how many planet pixels will fit in a landsat pixel
    row_fact = int(np.floor(planet_rows/oli_rows))
    col_fact = int(np.floor(planet_cols/oli_cols))
    
    # it is not necessarily an even number given the spatial extraction, so add in the modulus
    row_mod = planet_rows%oli_rows
    col_mod = planet_cols%oli_cols
    
    # get the ranges of pixels to grab from planet
    row_sets = np.arange(0,planet_rows, row_fact) + row_mod
    col_sets = np.arange(0, planet_cols, col_fact) + col_mod
    
    # sample via landsat
    samples = []
    for r in range(oli_rows-1):
        for c in range(oli_cols-1):
            
            # get the landsat pixel
            oli_pix_val_red = oli_imArr_red[r,c]
            oli_pix_val_nir = oli_imArr_nir[r,c]
            
            # get the pixels from planet if the pixel is not 0 or -9999
            if oli_pix_val_red is not np.ma.masked: 
                pl_row_start = row_sets[r]
                pl_row_end = row_sets[r+1]
                pl_col_start = col_sets[c]
                pl_col_end = col_sets[c+1]
                
                pl_pix_red = pl_imArr_red[pl_row_start:pl_row_end, pl_col_start:pl_col_end].ravel()
                pl_pix_nir = pl_imArr_nir[pl_row_start:pl_row_end, pl_col_start:pl_col_end].ravel()
                pl_ndvi = (pl_pix_nir - pl_pix_red) / (pl_pix_nir + pl_pix_red)
                pl_ndvi_mean = np.mean(pl_ndvi)
                pl_ndvi_std = np.std(pl_ndvi)
                pl_ndvi_max = np.max(pl_ndvi)
                pl_ndvi_min = np.min(pl_ndvi)
                
                oli_ndvi = (oli_pix_val_nir - oli_pix_val_red) / (oli_pix_val_nir + oli_pix_val_red)
                
                
                samples.append(np.array([[oli_ndvi], [pl_ndvi_mean], [pl_ndvi_std], [pl_ndvi_max - pl_ndvi_min]]))
    
    return samples