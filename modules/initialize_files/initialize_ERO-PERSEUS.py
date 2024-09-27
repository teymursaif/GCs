import os, sys
import numpy as np
import matplotlib.pyplot as plt

def initialize_params() :
    print (f"{bcolors.OKCYAN}- Initializing the pipeline ... "+ bcolors.ENDC)
    # Defining the working directories
    global working_directory,input_dir,output_dir,data_dir,main_data_dir,clean_data_dir,img_dir,sex_dir,fit_dir,plots_dir,\
    detection_dir,cats_dir,psfs_dir,art_dir,final_cats_dir,temp_dir,sbf_dir,psf_dir, check_plots_dir, external_dir, data_dir_orig, sub_data_dir
    # Getting the objcts and data info
    global gal_id, gal_name, ra, dec, distance, filters, comments, gal_params, gal_methods, gal_data_name
    # Configuring the pipeline parameters
    global PRIMARY_FRAME_SIZE_ARCSEC, FRAME_SIZE_ARCSEC, GAL_FRAME_SIZE_ARCSEC, N_ART_GCS, N_SIM_GCS, PSF_IMAGE_SIZE, INSTR_FOV, COSMIC_CLEAN, \
    PHOTOM_APERS, FWHMS_ARCSEC, APERTURE_SIZE, PSF_REF_RAD_FRAC, BACKGROUND_ANNULUS_START, BACKGROUND_ANNULUS_TICKNESS, TARGETS, APERTURE_SIZE, \
    MAG_LIMIT_CAT, CROSS_MATCH_RADIUS_ARCSEC, GC_SIZE_RANGE, GC_MAG_RANGE, RATIO_OVERSAMPLE_PSF, PSF_PIXEL_SCALE, PSF_SIZE, MODEL_PSF, \
    PIXEL_SCALES, ZPS, PRIMARY_FRAME_SIZE, FRAME_SIZE, GAL_FRAME_SIZE, EXPTIME, GAIN, GC_REF_MAG, PSF_PIXELSCL_KEY, FWHM_LIMIT, INPUT_ZP, INPUT_EXPTIME, INPUT_GAIN, \
    MAG_LIMIT_SAT, MAG_LIMIT_PSF, GC_SEL_PARAMS, ELL_LIMIT_PSF, GC_SIM_MODE, MERGE_CATS, MERGE_SIM_GC_CATS, MERGE_GC_CATS, EXTRACT_DWARFS,\
    PARAM_SEL_METHOD, PARAM_SEL_RANGE, EXTERNAL_CROSSMATCH, EXTERNAL_CROSSMATCH_CAT, PARAM_ADD_UPPER, PARAM_ADD_LOWER, APER_SIZE_IN_FWHM,\
    PARAM_PERCENTILE_LOWER, PARAM_PERCENTILE_UPPER, PARAM_DESC, MEDIAN_FILTER_SIZE
    global SE_executable,galfit_executable,swarp_executable

    FWHMS_ARCSEC = {}
    APERTURE_SIZE = {}
    PSF_REF_RAD_FRAC = {}
    PIXEL_SCALES = {}
    ZPS = {}
    PRIMARY_FRAME_SIZE = {}
    FRAME_SIZE = {}
    GAL_FRAME_SIZE = {}
    EXPTIME = {}
    GAIN = {}
    INPUT_ZP = {}
    INPUT_EXPTIME = {}

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################
    ##### PARAMETERS THAT USER NEEDS TO CONFIGURE

    #-------------------------- SETUP AND DATA PREPRATION --------------------------

    ### (if ZP, EXPTIME and GAIN are missing from the header, define them for a given filter)

    WORKING_DIR = './'
    main_data_dir = WORKING_DIR+'ERO-data/ERO-PERSEUS/'
    #PRIMARY_FRAME_SIZE_ARCSEC = 1200 #arcsec
    #FRAME_SIZE_ARCSEC = 1200 #cut-out size from the original frame for the general anlaysis (arcsec)

    # defining the executables (what you type in the command-line that executes the program)
    SE_executable = 'sex'
    swarp_executable = 'swarp'
    #SE_executable = 'sextractor'
    #swarp_executable = 'SWarp'

    ### (if ZP, EXPTIME and GAIN are missing from the header, define them for a given filter)
    INPUT_ZP = {'VIS-DET':30.1,'VIS':30.1,'NISP-Y':30,'NISP-J':30,'NISP-H':30}
    INPUT_EXPTIME = {'VIS-DET':1,'VIS':1,'NISP-Y':1,'NISP-J':1,'NISP-H':1}
    INPUT_GAIN = {'VIS-DET':2,'VIS':2,'NISP-Y':1,'NISP-J':1,'NISP-H':1}

    # ------------------------------ GALAXIES/TARGETS ------------------------------

    # List of targets as a string with:
    # Object-ID Object name RA Dec Distance-in-Mpc List-of-filters comment
    # About list of filters: first filter is the detection filter (separated by ",")
    # comments: LSB,N,etc
    # (lines with # in the beginning will be skipped)
    # example: '1 DF44 195.2416667 +26.9763889 100 F814W,F475W,F606W'

    # HST
    #TARGETS = ['1 MATLAS2019 226.33460 +01.81282 25 HST-ACS-F606W,HST-ACS-F814W LSB,nN']
    #GC_REF_MAG = {'HST-ACS-F606W':-7.5,'HST-ACS-F814W':-8.0}

    #JWST
    #TARGETS = ['1 CEERS-LSB1 214.8588333333 +52.7629166667 80 F115W,F150W,F200W,F277W,F356W,F444W LSB,SF']
    #GC_REF_MAG = {'F115W':-8.0,'F150W':-8.0,'F200W':-8.0,'F277W':-8.0,'F356W':-8.0,'F444':-8.0}

    #Euclid SIMS
    #TARGETS = ['1 EUC-SIM1 231.50075 +30.45227 20 VIS E,N']
    #TARGETS = ['2 DWARF-MER-SIM 269.06658 +65.00640 20 VIS LSB,N']
    #GC_REF_MAG = {'VIS':-8}

    #Euclid ERO
    TARGETS = []

    from astropy.wcs import WCS
    from astropy.io import fits
    frame = main_data_dir+'/ERO-PERSEUS_VIS.fits'
    main = fits.open(frame)
    data = main[0].data
    hdr = main[0].header
    X = hdr['NAXIS1']
    Y = hdr['NAXIS2']
    S = np.max([X,Y])
    W=WCS(frame)
    N = 20
    PRIMARY_FRAME_SIZE_ARCSEC = int(S/N/10)+5
    FRAME_SIZE_ARCSEC = int(S/N/10)+5
    FRAME_SIZE_PIXELS = FRAME_SIZE_ARCSEC*10
    print (FRAME_SIZE_ARCSEC)
    for j in range(N):
        for i in range(N):
            m = (j*N+i)
            xc =  (i) * (X/N) + (X/N)/2.
            yc =  (j) * (Y/N) + (Y/N)/2.
            ra, dec = W.all_pix2world(yc,xc,0)
            #if m in [3,4,5,7,0,2,6,8] : #[1,3,4,5,7]  [0,2,6,8] 
            if m>=-1: #0,1,2,3,4,5,6,7,8
                target_str = str(m)+' ERO-PERSEUS ERO-PERSEUS-'+str(m)+' '+str(ra)+' '+str(dec)+' 73 VIS,NISP-Y,NISP-J,NISP-H MAKE_GC_CAT ---' #MAKE_CAT,SIM_GC, ,NISP-Y,NISP-J,NISP-H
                
                print ('checking freame '+str(m)+'...')
                w = FRAME_SIZE_PIXELS/2-10
                cropped_data = data[int(xc-w):int(xc+w),int(yc-w):int(yc+w)]
                #plt.imshow(cropped_data,vmin=-10,vmax=10)
                #plt.savefig('test.png')
                cropped_data_nonzero = cropped_data[abs(cropped_data)>1e-10]
                N_nenzero_pixels = len(cropped_data_nonzero)
                #print (median_pixel_value)
                
                if (N_nenzero_pixels) > 25000 :
                    print (N_nenzero_pixels)
                    TARGETS.append([target_str])
                else:
                    print (N_nenzero_pixels)
                    print (f"{bcolors.FAIL}*** Frame looks empty, skipping the frame..."+ bcolors.ENDC)
            #print (target_str)

    print (TARGETS)

    MERGE_CATS = False
    MERGE_SIM_GC_CATS = False
    MERGE_GC_CATS = True

    # NOTE: possible methods -> RESAMPLE_DATA, MODEL_PSF, FIT_GAL, USE_SUB_GAL, MAKE_CAT, MAKE_GC_CAT
    # NOTE: possible comments -> MASSIVE,DWARF,LSB

    global TABLES
    TABLES = {}
    #TABLES['acsfcs']='./archival_tables/ACS-FCS-GCs.fits'
    #TABLES['fornax-spec-gcs']='./archival_tables/Fornax_spec_UCDs_and_GCs.fits'
    #TABLES['fornax-spec-stars']='./archival_tables/Fornax_spec_foreground_stars.fits'

    # ------------------------------  GALAXY FITTING ------------------------------

    GAL_FRAME_SIZE_ARCSEC  = 1*FRAME_SIZE_ARCSEC  #cut-out size from the original frame for sersic fitting anlaysis (arcsec)

    # ---------------------- SOURCE DETECTION AND PHOTOMETRY ----------------------

    MEDIAN_FILTER_SIZE = 3
    PHOTOM_APERS = '1,2,4,8,12,16,20,30,40' #aperture-sizes (diameters) in pixels for aperture photometry with Sextractor
    BACKGROUND_ANNULUS_START = 5 #The size of background annulus for forced photoemtry as a factor of FWHM
    BACKGROUND_ANNULUS_TICKNESS = 20 # the thickness of the background annulus in pixels
    CROSS_MATCH_RADIUS_ARCSEC = 0.25
    APER_SIZE_IN_FWHM = {'VIS-DET':1.5,'VIS':1.5,'NISP-Y':1.5,'NISP-J':1.5,'NISP-H':1.5}
    MAG_LIMIT_CAT = 28.0
    EXTRACT_DWARFS = False

    # -------------------------------- PSF MODELING -------------------------------

    ### if PSF is given by the user (in the "psf_dir" directory)
    PSF_PIXELSCL_KEY = 'PIXELSCL'
    PSF_PIXEL_SCALE = 0.0 #if 'PIXELSCL' is not in the header, specify it here.
    ### for making PSF (method=MODEL_PSF)
    MODEL_PSF = True
    RATIO_OVERSAMPLE_PSF = 10 #do not go beyond 10, this will have consequences for undersampling later
    PSF_IMAGE_SIZE = 100 #PSF size in the instruments pixel-scale
    MAG_LIMIT_PSF = 21
    MAG_LIMIT_SAT = 19
    ELL_LIMIT_PSF = 0.1
    #FWHM_UPPER_LIMIT_PSF =
    #FWHM_LOWER_LIMIT_PSF =


    #------------------------------ GC SIMULATION ------------------------------
    N_ART_GCS = 670
    N_SIM_GCS = 1
    COSMIC_CLEAN = False #does not work at the moment anyways...
    GC_SIZE_RANGE = [1,1.01] #lower value should be small enough to make some point-sources for performance check, in pc
    GC_MAG_RANGE = [-12,-5]
    GC_REF_MAG = {'VIS-DET':-8,'VIS':-8,'NISP-Y':-8.45,'NISP-J':-8.45,'NISP-H':-8.45 } #magnitude of a typical GC in the given filters should be defined here.
    GC_SIM_MODE = 'UNIFORM' # 'UNIFORM' or 'CONCENTRATED'

    #------------------------------ GC SELECTION -------------------------------

    GC_SEL_PARAMS = ['CI_2_4_VIS','ELLIPTICITY_VIS',['F_MAG_APER_CORR_VIS','F_MAG_APER_CORR_NISP-Y'],['F_MAG_APER_CORR_NISP-Y','F_MAG_APER_CORR_NISP-J'],\
                        ['F_MAG_APER_CORR_NISP-J','F_MAG_APER_CORR_NISP-H']]#,'CI_4_8','CI_8_12']#,'CI_2_4','CI_4_6','CI_6_8','CI_8_10','CI_10_12','ELLIPTICITY']
    PARAM_DESC = ['compactness','ellipticity','color','color','color']
    PARAM_ADD_LOWER = [-0.1,-0.1,-0.3,-0.15,-0.15]
    PARAM_ADD_UPPER = [-0.1,0.0,0.3,0.15,0.15]
    PARAM_PERCENTILE_LOWER = [5,1,1,1,1]
    PARAM_PERCENTILE_UPPER = [95,99,99,99,99]


    EXTERNAL_CROSSMATCH = True
    EXTERNAL_CROSSMATCH_CAT = './archival_tables/PERSEUS/ERO-PERSEUS-CFHT.fits'

    PARAM_SEL_METHOD = 'MANUAL'
    PARAM_SEL_RANGE = {'ELLIPTICITY':[-0.01,1.01],'F_MAG_APER_CORR_VIS':[22,28],'color0':['VIS','VIS',-0.5,0.5],'CI_2_4_VIS':[-10,10]}#,\
    #'color1':['VIS','NISP-Y',-2,3],'color2':['NISP-Y','NISP-J',-2,2],'color3':['NISP-J','NISP-H',-2,2]}#,\
   # 'color5':['u','i',1.0,4.0], 'color6':['g','i',0.5,2.5], 'color7':['r','i',0,1], 'color8':['i','k',1,3.5]}   # clean selection


    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    # Don't mind this part! (!!! DO NOT CHANGE UNLESS YOU KNOW WHATYOU ARE DOING)

    working_directory = WORKING_DIR

    input_dir = working_directory+'inputs/'
    output_dir = working_directory+'outputs/'

    data_dir = input_dir+'data/'
    psf_dir = input_dir+'psf/'
    clean_data_dir = output_dir+'clean_data/'
    sub_data_dir = output_dir+'sub_data/'
    img_dir = output_dir+'img/'
    sex_dir = output_dir+'sex/'
    fit_dir = output_dir+'fit/'
    plots_dir = output_dir+'plots/'
    detection_dir = output_dir+'detection/'
    cats_dir = output_dir+'cats/'
    psfs_dir = output_dir+'psfs/'
    art_dir = output_dir+'artificial/'
    final_cats_dir = output_dir+'final_cats/'
    temp_dir = output_dir+'temp_files/'
    sbf_dir = output_dir+'sbf/'
    check_plots_dir = output_dir+'check_plots/'
    external_dir = input_dir+'external/'
    data_dir_orig = data_dir
    galfit_executable = external_dir+'galfit'

    global log_file
    log_file = './perseus.LOG'

    for dir in [working_directory,input_dir,output_dir,data_dir,main_data_dir,clean_data_dir,img_dir,sex_dir,fit_dir,plots_dir,\
    detection_dir,cats_dir,psfs_dir,art_dir,final_cats_dir,temp_dir,sbf_dir,psf_dir,check_plots_dir,sub_data_dir] :
        if not os.path.exists(dir): os.makedirs(dir)

    gal_params = {}
    gal_methods = {}
    gal_data_name = {}
    for line in TARGETS:
        #print (line)
        line = line[0].split(' ')
        gal_id, data_name, gal_name, ra, dec, distance, filters, methods, comments = int(line[0]), str(line[1]), str(line[2]),\
        float(line[3]), float(line[4]), float(line[5]), np.array(line[6].split(',')), np.array(line[7].split(',')), np.array((line[8].split('\n')[0]).split(','))
        gal_params[gal_id] = [gal_name, ra, dec, distance, filters, comments]
        gal_methods[gal_id] = methods
        gal_data_name[gal_id] = data_name
        #print (gal_name, ra, dec, distance, filters, comments)


############################################################

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

def welcome():

    #print (f"\n{bcolors.OKCYAN}   *****************************************"+ bcolors.ENDC)
    print (f"{bcolors.OKCYAN} \n+ GCTOOLS (version: August 2023) "+ bcolors.ENDC)
    print (f"+ Developed by Teymoor Saifollahi "+ bcolors.ENDC)
    print (f"+ Kapteyn Astronomical Institute"+ bcolors.ENDC)
    print (f"+ contact: saifollahi@astro.rug.nl\n"+ bcolors.ENDC)
    #print (f"{bcolors.OKCYAN}   *****************************************\n"+ bcolors.ENDC)

############################################################

def finalize(gal_id):
    rm_keys = ['*.fits','*.log','galfit*','*.xml']
    for rm_key in rm_keys:
        os.system('rm '+rm_key)
    print (f"{bcolors.OKGREEN}- Everything is done for this objects."+ bcolors.ENDC)


############################################################

welcome()
initialize_params()

#export PATH="/net/cannon/data/users/saifollahi/miniconda3/bin:
#/net/cannon/data/users/saifollahi/miniconda3/condabin:
#/net/cannon/usr/lib64/qt-3.3/bin:/net/cannon/usr/local/bin:
#/net/cannon/usr/local/sbin:/net/cannon/usr/bin:/net/cannon/usr/sbin
#:/net/cannon/bin:/sbin:$PATH"
