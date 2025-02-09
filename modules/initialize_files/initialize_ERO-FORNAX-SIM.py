import os, sys
import numpy as np

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
    PIXEL_SCALES, ZPS, PRIMARY_FRAME_SIZE, FRAME_SIZE, GAL_FRAME_SIZE, EXPTIME, GAIN, GC_REF_MAG, PSF_PIXELSCL_KEY, FWHM_LIMIT, INPUT_ZP, INPUT_EXPTIME, \
    MAG_LIMIT_SAT, MAG_LIMIT_PSF, GC_SEL_PARAMS, ELL_LIMIT_PSF, GC_SIM_MODE, MERGE_CATS, MERGE_SIM_GC_CATS, MERGE_GC_CATS, EXTRACT_DWARFS,\
    PARAM_SEL_METHOD, PARAM_SEL_RANGE, EXTERNAL_CROSSMATCH, EXTERNAL_CROSSMATCH_CAT, MAG_LIMIT_FILTER 
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
    PRIMARY_FRAME_SIZE_ARCSEC = 1200 #arcsec
    FRAME_SIZE_ARCSEC = 1200 #cut-out size from the original frame for the general anlaysis (arcsec)

    # defining the executables (what you type in the command-line that executes the program)
    SE_executable = 'sex'
    swarp_executable = 'swarp'
    #SE_executable = 'sextractor'
    #swarp_executable = 'SWarp'

    ### (if ZP, EXPTIME and GAIN are missing from the header, define them for a given filter)
    INPUT_ZP = {'VIS':30,'NISP-Y':30,'NISP-J':30,'NISP-H':30}
    INPUT_EXPTIME = {'VIS':565,'NISP-Y':81,'NISP-J':81,'NISP-H':81}
    INPUT_GAIN = {'VIS':2,'NISP-Y':1,'NISP-J':1,'NISP-H':1}

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

    TARGETS.append(['1 ERO-FORNAX ERO-FORNAX-1 54.41498968710675 -35.58635104214481 20 VIS,NISP-Y,NISP-J,NISP-H MAKE_GC_CAT ---']) 
    TARGETS.append(['2 ERO-FORNAX ERO-FORNAX-2 54.020110517026325 -35.58700320641283 20 VIS,NISP-Y,NISP-J,NISP-H MAKE_GC_CAT ---'])
    TARGETS.append(['3 ERO-FORNAX ERO-FORNAX-3 53.625231108456134 -35.58636741821876 20 VIS,NISP-Y,NISP-J,NISP-H MAKE_GC_CAT ---'])
    TARGETS.append(['4 ERO-FORNAX ERO-FORNAX-4 54.41340800663024 -35.2638620852018 20 VIS,NISP-Y,NISP-J,NISP-H MAKE_GC_CAT ---'])

    TARGETS.append(['5 ERO-FORNAX ERO-FORNAX-5 54.02010052441082 -35.26450654419051 20 VIS,NISP-Y,NISP-J,NISP-H MAKE_GC_CAT ---'])
    TARGETS.append(['6 ERO-FORNAX ERO-FORNAX-6 53.62679280653799 -35.263878267794965 20 VIS,NISP-Y,NISP-J,NISP-H MAKE_GC_CAT ---'])
    TARGETS.append(['7 ERO-FORNAX ERO-FORNAX-7 54.41183879946199 -34.94135934996358 20 VIS,NISP-Y,NISP-J,NISP-H MAKE_GC_CAT ---'])
    TARGETS.append(['8 ERO-FORNAX ERO-FORNAX-8 54.020090610601954 -34.9419961237686 20 VIS,NISP-Y,NISP-J,NISP-H MAKE_GC_CAT ---'])

    TARGETS.append(['9 ERO-FORNAX ERO-FORNAX-9 53.62834218888022 -34.94137533958052 20 VIS,NISP-Y,NISP-J,NISP-H MAKE_GC_CAT ---'])

    MERGE_CATS = True
    MERGE_SIM_GC_CATS = True
    MERGE_GC_CATS = True

    # NOTE: possible methods -> RESAMPLE_DATA, MODEL_PSF, FIT_GAL, USE_SUB_GAL, MAKE_CAT, MAKE_GC_CAT
    # NOTE: possible comments -> MASSIVE,DWARF,LSB

    global TABLES
    TABLES = {}
    TABLES['acsfcs']='./archival_tables/ACS-FCS-GCs.fits'
    TABLES['fornax-spec-gcs']='./archival_tables/Fornax_spec_UCDs_and_GCs.fits' #only Saifollahi+2021b
    TABLES['fornax-spec-gcs-all']='./archival_tables/Fornax_spec_UCDs_and_GCs_all.fits'
    TABLES['fornax-spec-stars']='./archival_tables/Fornax_spec_foreground_stars.fits'
    TABLES['fornax-spec-galaxies']='./archival_tables/Fornax_spec_background_galaxies.fits'
    TABLES['gaia-stars']='./archival_tables/gaia_dr3_sources.fits'

    # ------------------------------  GALAXY FITTING ------------------------------

    GAL_FRAME_SIZE_ARCSEC  = 1*FRAME_SIZE_ARCSEC  #cut-out size from the original frame for sersic fitting anlaysis (arcsec)

    # ---------------------- SOURCE DETECTION AND PHOTOMETRY ----------------------

    PHOTOM_APERS = '1,2,4,8,12,16,20,30,40' #aperture-sizes (diameters) in pixels for aperture photometry with Sextractor
    BACKGROUND_ANNULUS_START = 4 #The size of background annulus for forced photoemtry as a factor of FWHM
    BACKGROUND_ANNULUS_TICKNESS = 20 # the thickness of the background annulus in pixels
    CROSS_MATCH_RADIUS_ARCSEC = 0.25
    MAG_LIMIT_CAT = 25.0 #25
    EXTRACT_DWARFS = False

    # -------------------------------- PSF MODELING -------------------------------

    ### if PSF is given by the user (in the "psf_dir" directory)
    PSF_PIXELSCL_KEY = 'PIXELSCL'
    PSF_PIXEL_SCALE = 0.0 #if 'PIXELSCL' is not in the header, specify it here.
    ### for making PSF (method=MODEL_PSF)
    MODEL_PSF = True
    RATIO_OVERSAMPLE_PSF = 10 #do not go beyond 10, this will have consequences for undersampling later
    PSF_IMAGE_SIZE = 40 #PSF size in the instruments pixel-scale
    MAG_LIMIT_PSF = 21 #19 for NISP
    MAG_LIMIT_SAT = 19 #17 for NISP #saturation limit
    ELL_LIMIT_PSF = 0.1
    #FWHM_UPPER_LIMIT_PSF =
    #FWHM_LOWER_LIMIT_PSF =


    #------------------------------ GC SIMULATION ------------------------------
    N_ART_GCS = 250
    N_SIM_GCS = 1
    COSMIC_CLEAN = False #does not work at the moment anyways...
    GC_SIZE_RANGE = [2,8] #lower value should be small enough to make some point-sources for performance check, in pc
    GC_MAG_RANGE = [-10,-5]
    GC_REF_MAG = {'VIS':-8, 'NISP-Y':-8.5,'NISP-J':-8.5,'NISP-H':-8.5 } #magnitude of a typical GC in the given filters should be defined here.
    GC_SIM_MODE = 'UNIFORM' # 'UNIFORM' or 'CONCENTRATED'

    #------------------------------ GC SELECTION -------------------------------

    GC_SEL_PARAMS = ['CI_2_4','CI_4_8','CI_8_12']#,'CI_2_4','CI_4_6','CI_6_8','CI_8_10','CI_10_12','ELLIPTICITY']
    EXTERNAL_CROSSMATCH = True
    EXTERNAL_CROSSMATCH_CAT = './archival_tables/ERO-FDS-ugriJKs.fits'

    PARAM_SEL_METHOD = 'MANUAL'
    PARAM_SEL_RANGE = {'color1':['VIS','NISP-Y',-0.1,0.7],'color2':['NISP-Y','NISP-J',-0.1,0.5],'color3':['NISP-J','NISP-H',-0.1,0.5], \
        'color5':['u','i',1.5,3.5], 'color6':['g','i',0.6,1.4], 'color7':['r','i',0,0.6], 'color8':['i','k',1,3.5], \
        'ELLIPTICITY':[0,0.5],'F_MAG_APER_CORR':[15,30]}

    #PARAM_SEL_RANGE = {'color1':['VIS','NISP-Y',-0.4,1.2],'color2':['NISP-Y','NISP-J',-0.5,0.5],'color3':['NISP-J','NISP-H',-0.3,0.4], \
    #    'color4':['u','NISP-Y',1,4.5], 'color5':['g','NISP-Y',0,2.5], 'color6':['r','NISP-Y',-0.5,1.5], 'color7':['i','NISP-Y',-0.7,1], \
    #    'color8':['u','i',1,4.0], 'color9':['i','k',0.5,3.5], \
    #    'ELLIPTICITY':[0,0.5]}

    ####################################################################################################
    ####################################################################################################
    ####################################################################################################

    # Don't mind this part! (!!! DO NOT CHANGE UNLESS YOU KNOW WHATYOU ARE DOING)

    working_directory = WORKING_DIR

    input_dir = working_directory+'inputs/'
    output_dir = working_directory+'outputs/'
    main_data_dir = working_directory+'ERO-data/ERO-FORNAX/'#input_dir+'main_data/'

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
