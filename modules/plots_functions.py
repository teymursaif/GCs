import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from modules.pipeline_functions import *
from modules.initialize import *
from modules.source_det import *

#plotting parameters
plt.style.use('tableau-colorblind10')
plt.rc('axes', labelsize=18)
plt.rc('xtick', labelsize=15)
plt.rc('ytick', labelsize=15)
plt.rc('axes', linewidth=1.7)
plt.rc('font',**{'family':'serif','serif':['Times']})
plt.rcParams['hatch.linewidth'] = 2.0
plt.rc('text', usetex=True)

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

############################################################

def assess_GC_simulations(gal_id):
    data_name = gal_data_name[gal_id]
    print (f"{bcolors.OKCYAN}- Making assessment plots for source detection and gc simulations"+ bcolors.ENDC)

    # for fornax ERO
    assess_GC_simulations_general(gal_id)
    #assess_GC_simulations_compactness(gal_id)

############################################################

def shorten_filter_name(filtername):
    short_filternames = ['F606W','F814W','VIS','F475W','NISP-Y','NISP-J','NISP-H','DES-g','DES-r','DES-i']
    for short_filtername in short_filternames:
        if short_filtername in filtername :
            return short_filtername

############################################################


def assess_GC_simulations_general(gal_id):

    gal_name, ra, dec, distance, filters, comments = gal_params[gal_id]
    data_name = gal_data_name[gal_id]

    # load main cats for plotting
    fn_det = filters[0]
    sfn_det = shorten_filter_name(fn_det)
    ART_GCs_cat = art_dir+gal_name+'_'+fn_det+'_ALL_ART_GCs.fits'
    DET_ART_GCs_cat = art_dir+gal_name+'_'+fn_det+'_ALL_DET_ART_GCs.fits'
    source_cat = cats_dir+gal_name+'_master_cat_forced.fits'
    art_gcs = (fits.open(ART_GCs_cat))[1].data
    det_art_gcs = (fits.open(DET_ART_GCs_cat))[1].data
    sources = (fits.open(source_cat))[1].data

    colors = ['red','green','cyan','blue','purple','black']

    if len(filters)==1:
        fig, ax = plt.subplots(6, len(filters)+1, figsize=(6*(len(filters)+1),7*3), constrained_layout=True)
    else:
        fig, ax = plt.subplots(6, len(filters), figsize=(6*len(filters),7*3), constrained_layout=True)

    n_filter = -1
    for fn in filters:
        try:
            sfn = shorten_filter_name(fn)
        except:
            sfn = fn
        n_filter = n_filter+1
        #1. mag differance (dmag - mag) in all filters
        mag_source_obs = det_art_gcs['F_MAG_APER_CORR_'+fn]
        mag_source_sim = det_art_gcs['GC_MAG_'+fn]
        ax[0][n_filter].scatter(mag_source_sim,mag_source_obs- mag_source_sim,s=20,color=colors[n_filter],alpha=0.5,marker='o',label='Detected Artificial GCs')
        #ax[0][n_filter].set_title(gal_name+', '+sfn)
        ax[0][n_filter].set_xlabel('$m_{'+sfn+', SIM} \ [mag]$')
        ax[0][n_filter].set_ylabel('$m_{'+sfn+', AP} - m_{'+sfn+', SIM} \ [mag]$')
        mag1 = GC_MAG_RANGE[0]+5*np.log10(distance*1e+5)-0.25+color(fn,fn_det)
        #mag2 = GC_MAG_RANGE[1]+5*np.log10(distance*1e+5)+0.25+color(fn,fn_det)
        mag2 =np.nanmax(mag_source_obs)
        ax[0][n_filter].set_xlim([mag1,mag2])
        ax[0][n_filter].set_ylim([-0.5,0.5])
        ax[0][n_filter].legend(loc='upper left')
        ax[0][n_filter].tick_params(which='both',direction='in')

        #2. color differance between filter 1 and 2 (if more than 1 filter is there)
        mag_source_obs = det_art_gcs['F_MAG_APER_CORR_'+fn]
        mag_det_source_obs = det_art_gcs['MAG_APER_CORR_'+fn_det]
        ax[1][n_filter].scatter(mag_source_obs,mag_source_obs-mag_det_source_obs,s=20,color=colors[n_filter],alpha=0.5,marker='o',label='Detected Artificial GCs')
        #ax[1][n_filter].set_title(gal_name+', '+sfn)
        ax[1][n_filter].set_xlabel('$m_{'+sfn+', SE-AP} \ [mag]$')
        ax[1][n_filter].set_ylabel('$m_{'+sfn+', SE-AP} - m_{'+sfn_det+', FP-AP} \ [mag]$')
        mag1 = GC_MAG_RANGE[0]+5*np.log10(distance*1e+5)-0.25+color(fn,fn_det)
        #mag2 = GC_MAG_RANGE[1]+5*np.log10(distance*1e+5)+0.25+color(fn,fn_det)
        mag2 =np.nanmax(mag_source_obs)
        ax[1][n_filter].set_xlim([mag1,mag2])
        ax[1][n_filter].set_ylim([color(fn,fn_det)-0.2,color(fn,fn_det)+0.2])
        ax[1][n_filter].legend(loc='upper left')
        ax[1][n_filter].tick_params(which='both',direction='in')

        #3. 4. fwhm and CI of the simulated GCs vs other objects in the real data (possibility to add external cats? later!)
        if fn == fn_det:
            param1 = det_art_gcs['FWHM_IMAGE_'+fn]
            param2 = det_art_gcs['MAG_AUTO_'+fn]
            param3 = sources['FWHM_IMAGE_'+fn]
            param4 = sources['MAG_AUTO_'+fn]
            ax[2][0].scatter(param4,param3,s=20,color='grey',alpha=0.2,marker='o',label='Detected Sources')
            ax[2][0].scatter(param2,param1,s=50,color=colors[n_filter],alpha=0.2,marker='o',label='Detected Artificial GCs')
            #ax[2][0].set_title(gal_name+', '+sfn)
            ax[2][0].set_xlabel('$m_{'+sfn+', AUTO} \ [mag]$')
            ax[2][0].set_ylabel('$FWHM_{'+sfn+', AUTO}\  [pixels]$')
            mag1 = GC_MAG_RANGE[0]+5*np.log10(distance*1e+5)-0.25+color(fn,fn_det)
            #mag2 = GC_MAG_RANGE[1]+5*np.log10(distance*1e+5)+0.25+color(fn,fn_det)
            mag2 = np.nanmax(param2)
            ax[2][0].set_xlim([mag1,mag2])
            fwhm1 = FWHMS_ARCSEC[fn]/PIXEL_SCALES[fn]*0.5
            fwhm2 = FWHMS_ARCSEC[fn]/PIXEL_SCALES[fn]*2.5
            ax[2][0].set_ylim([fwhm1,fwhm2])
            #ax[2][0].set_ylim([0,10])
            ax[2][0].legend(loc='upper left')
            ax[2][0].tick_params(which='both',direction='in')


        if fn == fn_det:
            for i in range(len(GC_SEL_PARAMS)):

                if PARAM_DESC[i] == 'color' :
                    continue
                    
                if i == len(filters):
                    break

                param1 = det_art_gcs[GC_SEL_PARAMS[i]] #+'_'+fn
                param2 = det_art_gcs['MAG_AUTO_'+fn]
                param3 = sources[GC_SEL_PARAMS[i]] #+'_'+fn
                param4 = sources['F_MAG_APER_CORR_'+fn]
                ax[3][i].scatter(param3,param4,s=20,color='grey',alpha=0.2,marker='o',label='Detected Sources')
                ax[3][i].scatter(param1,param2,s=50,color=colors[n_filter],alpha=0.2,marker='o',label='Detected Artificial GCs')
                #ax[2][n_filter+1].set_title(gal_name+', '+sfn)
                ax[3][i].set_ylabel('$m_{'+sfn+', AUTO} \ [mag]$')
                ax[3][i].set_xlabel('$'+(GC_SEL_PARAMS[i]).replace("_", "-")+'_{'+sfn+'} \ [mag]$')
                mag1 = GC_MAG_RANGE[0]+5*np.log10(distance*1e+5)-0.25+color(fn,fn_det)
                #mag2 = GC_MAG_RANGE[1]+5*np.log10(distance*1e+5)+0.25+color(fn,fn_det)
                mag2 = np.nanmax(param2)
                ax[3][i].set_ylim([mag2,mag1])
                ax[3][i].set_xlim([-0.2,1.5])
                ax[3][i].legend(loc='upper left')
                ax[3][i].tick_params(which='both',direction='in')

            ### add data from some tables

        #5. completness of detection in mag in DET filter
        #if fn == fn_det:
        mask = (art_gcs['ART_GC_FLAG_'+fn]>0)
        art_gcs = art_gcs[mask]
        mag_source_sim = art_gcs['GC_MAG_'+fn]

        mask = ( (det_art_gcs['F_MAG_APER_CORR_'+fn]>0) & (det_art_gcs['F_MAG_APER_CORR_'+fn]< 50) )
        det_art_gcs = det_art_gcs[mask]
        mag_source_obs = det_art_gcs['GC_MAG_'+fn]

        ax[4][n_filter].hist(mag_source_sim,bins=np.arange(10,30,0.2),color='black',alpha=1,label='ALL simulated GCs')
        ax[4][n_filter].hist(mag_source_obs,bins=np.arange(10,30,0.2),color=colors[n_filter],alpha=1,label='Detected simulated GCs')
        mag1 = GC_MAG_RANGE[0]+5*np.log10(distance*1e+5)-0.25+color(fn,fn_det)
        mag2 = GC_MAG_RANGE[1]+5*np.log10(distance*1e+5)+0.25+color(fn,fn_det)
        ax[4][n_filter].set_xlim([mag1,mag2])
        ax[4][n_filter].set_xlabel('$m_{'+sfn+', SIM} \ [mag]$')
        ax[4][n_filter].set_ylabel('N')
        ax[4][n_filter].legend(loc='upper left')

        #6. completness of detection in radial distance in all filters and the combination
        #if fn == fn_det:
        #r_source_sim = np.sqrt((art_gcs['RA_'+fn]-ra)**2 + (art_gcs['DEC_'+fn]-dec)**2)
        #r_source_obs = np.sqrt((det_art_gcs['RA_'+fn]-ra)**2 + (det_art_gcs['DEC_'+fn]-dec)**2)

        mask = (art_gcs['ART_GC_FLAG_'+fn]>0)
        art_gcs = art_gcs[mask]
        mag_source_sim = art_gcs['GC_MAG_'+fn]

        mask = ( (det_art_gcs['F_MAG_APER_CORR_'+fn]>0) & (det_art_gcs['F_MAG_APER_CORR_'+fn]< 50) )
        det_art_gcs = det_art_gcs[mask]
        mag_source_obs = det_art_gcs['GC_MAG_'+fn]

        r_source_sim = np.sqrt((art_gcs['RA_GC']-ra)**2 + (art_gcs['DEC_GC']-dec)**2)*3600
        r_source_obs = np.sqrt((det_art_gcs['RA_GC']-ra)**2 + (det_art_gcs['DEC_GC']-dec)**2)*3600

        r_source_sim = r_source_sim[mag_source_sim<(np.nanmax(mag_source_obs))]
        r_source_obs = r_source_obs[mag_source_obs<(np.nanmax(mag_source_obs))]

        ax[5][n_filter].hist(r_source_sim,bins=np.arange(0,np.nanmax(r_source_sim),np.nanmax(r_source_sim)/20),color='black',alpha=1,label='ALL simulated GCs')
        ax[5][n_filter].hist(r_source_obs,bins=np.arange(0,np.nanmax(r_source_sim),np.nanmax(r_source_sim)/20),color=colors[n_filter],alpha=1,label='Detected simulated GCs')
        #mag1 = GC_MAG_RANGE[0]+5*np.log10(distance*1e+5)-0.25+color(fn,fn_det)
        #mag2 = GC_MAG_RANGE[1]+5*np.log10(distance*1e+5)+0.25+color(fn,fn_det)
        #ax[3][i].set_xlim([mag1,mag2])
        ax[5][n_filter].set_xlabel('$R_{GAL, SIM} \ [arcsec]$')
        ax[5][n_filter].set_ylabel('N')
        ax[5][n_filter].legend(loc='upper left')

    fig.tight_layout()
    plt.savefig(plots_dir+gal_name+'_'+fn+'_check_plot.png',dpi=100)
    plt.close()

############################################################

def assess_GC_simulations_compactness(gal_id):

    gal_name, ra, dec, distance, filters, comments = gal_params[gal_id]
    data_name = gal_data_name[gal_id]

    # load main cats for plotting
    fn_det = filters[0]
    sfn_det = shorten_filter_name(fn_det)
    ART_GCs_cat = art_dir+gal_name+'_'+fn_det+'_ALL_ART_GCs.fits'
    DET_ART_GCs_cat = art_dir+gal_name+'_'+fn_det+'_ALL_DET_ART_GCs.fits'
    source_cat = cats_dir+gal_name+'_master_cat_forced.fits'
    art_gcs = (fits.open(ART_GCs_cat))[1].data
    det_art_gcs = (fits.open(DET_ART_GCs_cat))[1].data
    sources = (fits.open(source_cat))[1].data

    acsfcs_table = source_cat+'.+ACSFCS.fits'
    crossmatch(source_cat,TABLES['acsfcs'],'RA','DEC','RAJ2000','DEJ2000',2.*PIXEL_SCALES[fn_det],'',acsfcs_table)

    try :
        cat = acsfcs_table
        make_cat_topcat_friendly(cat,cat[:len(cat)-5]+'+.fits')
    except :
        donothing = 1

    acsfcs = (fits.open(acsfcs_table))[1].data
    acsfcs = acsfcs[acsfcs['pGC']>0.95]

    specgcs_table = source_cat+'.+SPECGCs.fits'
    crossmatch(source_cat,TABLES['fornax-spec-gcs'],'RA','DEC','RA','DEC',2.*PIXEL_SCALES[fn_det],'',specgcs_table)

    try :
        cat = specgcs_table
        make_cat_topcat_friendly(cat,cat[:len(cat)-5]+'+.fits')
    except :
        donothing = 1

    specgcs = (fits.open(specgcs_table))[1].data

    PHOTOM_APERS_ = PHOTOM_APERS.split(',')

    for i in range(len(PHOTOM_APERS_)-1):
        fig, ax = plt.subplots(figsize=(4,6), constrained_layout=True)
        aper1 = PHOTOM_APERS_[i]
        aper2 = PHOTOM_APERS_[i+1]
        param = 'CI_'+str(aper1)+'_'+str(aper2)+'_'+fn_det
        mag_param = 'MAG_APER_CORR_'+fn_det

        #ci_art_gcs = art_gcs[param]
        ci_det_art_gcs = det_art_gcs[param]
        ci_sources = sources[param]
        ci_acsfcs = acsfcs[param]
        ci_specgcs = specgcs[param]

        #mag_art_gcs = art_gcs[mag_param]
        mag_det_art_gcs = det_art_gcs[mag_param]
        mag_sources = sources[mag_param]
        mag_acsfcs = acsfcs[mag_param]
        mag_specgcs = specgcs[mag_param]

        ax.scatter(ci_sources,mag_sources,s=10,color='grey',alpha=0.1,marker='o',label='Detected Sources')
        ax.scatter(ci_acsfcs,mag_acsfcs,s=20,color='gold',alpha=1,marker='o',label='ACSFCS UCDs/GCs (pGC $>$ 0.95')
        ax.scatter(ci_det_art_gcs,mag_det_art_gcs,s=20,color='green',alpha=1,marker='o',label='Simulated GCs')
        ax.scatter(ci_specgcs,mag_specgcs,s=20,color='red',alpha=1,marker='o',label='Spectroscopic UCDs/GCs')

        ax.set_xlabel('$CI_{'+str(aper1)+'-'+str(aper2)+','+fn_det+'}$')
        ax.set_ylabel('$m_{'+str(fn_det)+'}$')

        ci1 = np.nanmedian(ci_det_art_gcs)-0.5
        ci2 = np.nanmedian(ci_det_art_gcs)+0.5

        ax.set_xlim([ci1,ci2])
        ax.set_ylim([MAG_LIMIT_CAT,MAG_LIMIT_SAT])

        plt.legend()
        plt.savefig(plots_dir+gal_name+'_'+fn_det+'_compactness-index-'+str(aper1)+'-'+str(aper2)+'.png',dpi=150)
        plt.close()

############################################################
