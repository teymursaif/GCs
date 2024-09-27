import os, sys
import math
import numpy as np
import matplotlib.pyplot as plt
#from astroquery.mast import Observations
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.ndimage import gaussian_filter
from scipy.ndimage import median_filter
from astropy.stats import sigma_clip
from astropy.visualization import *
from astropy.visualization import make_lupton_rgb
from astropy.table import Table, join_skycoord
from astropy import table
from astropy.stats import sigma_clip
import photutils
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.centroids import centroid_1dg, centroid_2dg, centroid_com, centroid_quadratic
#from photutils.utils import CutoutImage
import time as TIME
from modules.initialize import *
from modules.pipeline_functions import *
from scipy.ndimage import median_filter

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def select_gcs_for_param(param_sources,mag_sources,param_det_art_gcs,mag_det_art_gcs,add_upper,add_lower,percentile_upper,percentile_lower,keep_bad_values=0,label=''):
    mask = []

    #cat1 = open(cats_dir+'compactness_range_lower_limit_'+label+'.csv', 'w')
    #cat2 = open(cats_dir+'compactness_range_upper_limit_'+label+'.csv', 'w')

    mag_range = np.arange(18,28,0.1)
    upper_limits = np.zeros(len(mag_range))-99
    lower_limits = np.zeros(len(mag_range))-99

    #mask_obs = ((abs(mag_sources)<30) & (abs(param_sources)<30))
    mask_art = ((abs(mag_det_art_gcs)<30) & (abs(param_det_art_gcs)<30)) 
    
    #param_sources = param_sources[mask_obs]
    #mag_sources = mag_sources[mask_obs]

    param_det_art_gcs = param_det_art_gcs[mask_art]
    mag_det_art_gcs = mag_det_art_gcs[mask_art]

    plt.plot(param_det_art_gcs,mag_det_art_gcs,'k.',alpha=0.01)

    N = len(param_sources)

    for mag in mag_range:

        j = int((mag-np.min(mag_range))/0.1+0.5)

        #print (mag,j)

        m = mag_det_art_gcs[abs(mag_det_art_gcs-mag)<0.5]
        p = param_det_art_gcs[abs(mag_det_art_gcs-mag)<0.5]

        p = sigma_clip(p, 3, masked=False)

        if len(p) < 20:
            lower_limits[j] = -99
            upper_limits[j] = -99
        else:
            try : 
                lower_limits[j] = (np.percentile(p, percentile_lower))+add_lower
            except :
                lower_limits[j] = (np.percentile(p, 1))-0.0
            try :
            	upper_limits[j] = (np.percentile(p, percentile_upper))+add_upper
            except:
                upper_limits[j] = (np.percentile(p, 99.0))+0.0
    
    for mag in mag_range:

        j = int((mag-np.min(mag_range))/0.1+0.5)

        if lower_limits[j] < -90 :
            lower_limits[j] = np.nanmax(lower_limits)

        if upper_limits[j] < -90 :
            upper_limits[j] = np.nanmin(upper_limits[upper_limits>-90])

    
    lower_limits = median_filter(lower_limits, size=5, cval=0, mode='constant')
    upper_limits = median_filter(upper_limits, size=5, cval=0, mode='constant')

    plt.plot(lower_limits,mag_range,'k.')
    plt.plot(upper_limits,mag_range,'k.')
    #plt.xlim([-0.5,1.5])

    print_str = 'mag_range=['
    for iii in range(len(mag_range)):
        print_str = print_str + str(mag_range[iii])+','
    print_str=print_str[:len(print_str)-1]
    print_str = print_str +']'
    print (print_str)

    print_str = 'upper_limits=['
    for iii in range(len(mag_range)):
        print_str = print_str + str(upper_limits[iii])+','
    print_str=print_str[:len(print_str)-1]
    print_str = print_str +']'
    print (print_str)

    print_str = 'lower_limits=['
    for iii in range(len(mag_range)):
        print_str = print_str + str(lower_limits[iii])+','
    print_str=print_str[:len(print_str)-1]
    print_str = print_str +']'
    print (print_str)


    for i in range(N):

        mag = mag_sources[i]
        param = param_sources[i]

        j = int((mag-np.min(mag_range))/0.1+0.5)

        if (j<0) or (j>=len(mag_range)):
            mask.append(0)
            continue

        #if (mag>MAG_LIMIT_CAT) or (mag<MAG_LIMIT_SAT):
        #    mask.append(0)
        #    continue

        #mag0 = mag #int(mag*2+0.5)/2.

        #mag_range_mask = ((abs(mag_det_art_gcs-mag0)<2.0))
        #param_det_art_gcs_in_mag_range = param_det_art_gcs[mag_range_mask]

        param_lower_limit = lower_limits[j] #param_median-2.0*param_std-0.05
        param_upper_limit = upper_limits[j] #param_median+2.0*param_std+0.05

        if keep_bad_values == 0 :

            if (param > param_lower_limit) and (param < param_upper_limit) and (abs(mag) < 50):
                mask.append(1)
                plt.plot(param,mag,'g.',alpha=0.2)
            else :
                mask.append(0) #0
                plt.plot(param,mag,'r.',alpha=0.2)

        if keep_bad_values == 1 :

            if (param > param_lower_limit) and (param < param_upper_limit) and (abs(mag) < 50):
                mask.append(1)
                plt.plot(param,mag,'g.',alpha=0.2)

            elif (abs(param) > 50) and (abs(mag) < 50):
                mask.append(1)
                plt.plot(param,mag,'b.',alpha=0.2)

            else :
                mask.append(0) #0
                plt.plot(param,mag,'r.',alpha=0.2)



    plt.savefig(check_plots_dir+'parameter_range_'+label+'.png')
    plt.xlim([-1,1])
    plt.savefig(check_plots_dir+'parameter_range_'+label+'+xlim.png')
    plt.close()

    #cat1.close()
    #cat2.close()

    mask = np.array(mask)
    return mask

def select_GC_candidadates(gal_id):
    gal_name, ra, dec, distance, filters, comments = gal_params[gal_id]
    data_name = gal_data_name[gal_id]

    fn_det = filters[0]
    #ART_GCs_cat = art_dir+gal_name+'_'+fn_det+'_ALL_ART_GCs.fits'
    
    try :
        DET_ART_GCs_cat = art_dir+gal_name+'_'+fn_det+'_ALL_DET_ART_GCs.fits_' #_random
        det_art_gcs = (fits.open(DET_ART_GCs_cat))[1].data
    except:
        print ('*** GC simulations are not found. Pipeline try to use some older simulated data that seems relevant, if available ... ')
        DET_ART_GCs_cat = cats_dir+data_name+'_ALL_DET_ART_GCs_merged.fits'
        det_art_gcs = (fits.open(DET_ART_GCs_cat))[1].data


    source_cat = cats_dir+gal_name+'_master_cat_forced.fits'
    selected_gcs_cat = final_cats_dir+gal_name+'_selected_GCs.fits'
    shutil.copy(source_cat,selected_gcs_cat)

    #art_gcs = (fits.open(ART_GCs_cat))[1].data
    sources = (fits.open(source_cat))[1].data

    if 'DET' in fn_det:
        mag_param = 'F_MAG_APER_CORR_'+filters[1]
    else:
        mag_param = 'F_MAG_APER_CORR_'+fn_det
        
    mag_mask = (sources[mag_param] < GC_MAG_RANGE[1]+5*np.log10(distance*1e+5))
    sources = sources[mag_mask]

    selected_gcs_mask = np.ones(len(sources))

    N_param = len(GC_SEL_PARAMS)

    for i in range(N_param):
        param = GC_SEL_PARAMS[i]
        param_desc = PARAM_DESC[i]
        print (param)

        if param_desc == 'color' :

            param1 = param[0]
            param2 = param[1]
            param = param1+'-'+param2

            param1_det_art_gcs = det_art_gcs[param1]
            param1_sources = sources[param1]

            param2_det_art_gcs = det_art_gcs[param2]
            param2_sources = sources[param2]

            param_det_art_gcs = param1_det_art_gcs - param2_det_art_gcs
            param_sources = param1_sources - param2_sources

            mag_det_art_gcs = det_art_gcs[mag_param]
            mag_sources = sources[mag_param]

            add_upper = PARAM_ADD_UPPER[i]
            add_lower = PARAM_ADD_LOWER[i]
            percentile_upper = PARAM_PERCENTILE_UPPER[i]
            percentile_lower = PARAM_PERCENTILE_LOWER[i]
            mask = select_gcs_for_param(param_sources,mag_sources,param_det_art_gcs,mag_det_art_gcs,add_upper,add_lower,percentile_upper,percentile_lower,\
                   keep_bad_values=1,label=param+'+')
            selected_gcs_mask = selected_gcs_mask * mask

        else : 
            #param = param + '_' + fn_det
            #param_art_gcs = art_gcs[param]
            param_det_art_gcs = det_art_gcs[param]
            param_sources = sources[param]
            #mag_art_gcs = art_gcs[mag_param]
            mag_det_art_gcs = det_art_gcs[mag_param]
            mag_sources = sources[mag_param]
            add_upper = PARAM_ADD_UPPER[i]
            add_lower = PARAM_ADD_LOWER[i]
            percentile_upper = PARAM_PERCENTILE_UPPER[i]
            percentile_lower = PARAM_PERCENTILE_LOWER[i]
            mask = select_gcs_for_param(param_sources,mag_sources,param_det_art_gcs,mag_det_art_gcs,add_upper,add_lower,percentile_upper,percentile_lower,\
                   keep_bad_values=0,label=param+'+')
            selected_gcs_mask = selected_gcs_mask * mask

    selected_gcs_mask = selected_gcs_mask.astype(np.bool_)
    selected_gcs_data = sources[selected_gcs_mask]

    selected_gcs = fits.open(selected_gcs_cat)
    selected_gcs[1].data = selected_gcs_data
    selected_gcs.writeto(selected_gcs_cat, overwrite=True)

    cat = selected_gcs_cat
    make_cat_topcat_friendly(cat,cat[:len(cat)-5]+'+.fits')
    shutil.copy(cat[:len(cat)-5]+'+.fits',cat[:len(cat)-5]+'+.size-selected.fits')

    print ('--- number of selected GCs after filter(s) in compactness is: ', len(selected_gcs_data) )

    ### extra selection

    if PARAM_SEL_METHOD == 'MANUAL':

        if EXTERNAL_CROSSMATCH == True:
            selected_gcs_cat_cm = final_cats_dir+gal_name+'_selected_GCs.EXT-MATCHED.fits'
            crossmatch_left_wing(selected_gcs_cat,EXTERNAL_CROSSMATCH_CAT,'RA','DEC','RA','DEC',5.*PIXEL_SCALES[fn_det],fn_det,selected_gcs_cat_cm)
            source_cat = selected_gcs_cat_cm

        elif EXTERNAL_CROSSMATCH == False:
            source_cat = selected_gcs_cat

        sources = (fits.open(source_cat))[1].data

        for param in PARAM_SEL_RANGE.keys():
            #print (param)
            if 'color' in str(param):
                f1 = (PARAM_SEL_RANGE[param])[0]
                f2 = (PARAM_SEL_RANGE[param])[1]

                try:
                    param_f1 = 'F_MAG_APER_CORR_'+f1
                    param_f2 = 'F_MAG_APER_CORR_'+f2
                    param_min = (PARAM_SEL_RANGE[param])[2]
                    param_max = (PARAM_SEL_RANGE[param])[3]

                    try : m1 = sources[param_f1]
                    except: m1 = sources[f1]
                    try : m2 = sources[param_f2]
                    except: m2 = sources[f2]

                except:
                    param_f1 = 'MAG_APER_CORR_'+f1
                    param_f2 = 'F_MAG_APER_CORR_'+f2
                    param_min = (PARAM_SEL_RANGE[param])[2]
                    param_max = (PARAM_SEL_RANGE[param])[3]

                    try : m1 = sources[param_f1]
                    except: m1 = sources[f1]
                    try : m2 = sources[param_f2]
                    except: m2 = sources[f2]

                color = m1 - m2

                for i in range(len(color)) :
                    if (m1[i] > 0) and (m2[i] > 0) : 
                        donothing = 1
                    else: 
                        #print (m1[i], m2[i])
                        color[i] = 0.5*(param_min+param_max) 

                param_mask = ((color>=param_min) & (color<=param_max))
                sources = sources[param_mask]
                print ('--- number of selected GCs after filter in color ', f1, ' - ', f2, ' is: ', len(sources) )
                #expand_fits_table(source_cat,f1+'_'+f2,np.array(color))

            else :
                param_min = (PARAM_SEL_RANGE[param])[0]
                param_max = (PARAM_SEL_RANGE[param])[1]

                try:
                    param = param
                    param_mask = ((sources[param]>=param_min) & (sources[param]<=param_max))
                    sources = sources[param_mask]
                    print ('--- number of selected GCs after filter in ', param, ' is: ', len(sources))
                    
                except:
                    
                    param = param + '_' + fn_det
                    param_mask = ((sources[param]>=param_min) & (sources[param]<=param_max))
                    sources = sources[param_mask]
                    print ('--- number of selected GCs after filter in ', param, ' is: ', len(sources))

        selected_gcs_data = sources
        selected_gcs = fits.open(selected_gcs_cat)
        selected_gcs[1].data = sources
        selected_gcs.writeto(selected_gcs_cat, overwrite=True)

        cat = selected_gcs_cat
        make_cat_topcat_friendly(cat,cat[:len(cat)-5]+'+.fits')


