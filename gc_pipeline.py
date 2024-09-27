import os, sys
import warnings

copy_initialize_file_command = 'cp '+str(sys.argv[1])+' ./modules/initialize.py'
os.system(copy_initialize_file_command)
#print (sys.argv)

from modules.initialize import *
from modules.pipeline_functions import *
from modules.plots_functions import *
from modules.source_det import *
from modules.fit_galaxy import *
from modules.psf import *
from modules.gc_det import *
from datetime import datetime

os.system('rm '+log_file)

for gal_id in gal_params.keys():

    log = open(log_file,'a')

    gal_name, ra, dec, distance, filters, comments = gal_params[gal_id]
    methods = gal_methods[gal_id]

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    log.write(gal_name+'\n')
    log.write('START: '+dt_string+'\n')

    # step 0. inistialize pipleine and prepare data
    intro(gal_id)
    if 'PREPARE_DATA' in sys.argv:
        copy_data(gal_id)

    get_data_info(gal_id)

    if 'RESAMPLE'in methods:
        resample_data(gal_id)
    
    if 'PREPARE_FRAME' in sys.argv:
        make_galaxy_frames(gal_id)

    if 'MODEL_PSF' in methods:
        make_psf_all_filters(gal_id)

    initial_psf(gal_id)

    now0 = datetime.now()
    dt_string = now0.strftime("%d/%m/%Y %H:%M:%S")
    dt = (now0-now)
    log.write('- STEP 0 DONE: '+dt_string+', TIME SPENT: '+str(dt.seconds)+'\n')

    # step 1. fit sersic and subtract the light
    if 'FIT_GAL' in methods:
        fit_galaxy_sersic_all_filters(gal_id)

    now1 = datetime.now()
    dt = (now1-now0)
    dt_string = now1.strftime("%d/%m/%Y %H:%M:%S")
    log.write('- STEP 1 DONE: '+dt_string+', TIME SPENT: '+str(dt.seconds)+'\n')

    # step 2. detection, photometry and make the main source catalogue
    if 'MAKE_CAT' in methods:
        #try :
        make_source_cat_full(gal_id)
        #assess_photom(gal_id)
        #except:
        #log.write('*** Something gone wrong with ctalogue production, skipping the frame...')
        #print (f"{bcolors.FAIL}*** Something gone wrong with ctalogue production, skipping the frame (maybe frame is blank?)"+ bcolors.ENDC)

    now2 = datetime.now()
    dt = (now2-now1)
    dt_string = now2.strftime("%d/%m/%Y %H:%M:%S")
    log.write('- STEP 2 DONE: '+dt_string+', TIME SPENT: '+str(dt.seconds)+'\n')

    # step 3. GC analysis: completeness, selection, measurments
    if 'SIM_GC' in methods:
        simulate_GCs_all(gal_id)
        make_source_cat_for_sim(gal_id)
        #assess_GC_simulations(gal_id)

    now3 = datetime.now()
    dt = (now3-now2)
    dt_string = now3.strftime("%d/%m/%Y %H:%M:%S")
    log.write('- STEP 3 DONE: '+dt_string+', TIME SPENT: '+str(dt.seconds)+'\n')

    # step 4. GC selection, assessment, GC catalogs and properties
    if 'MAKE_GC_CAT' in methods:
        select_GC_candidadates(gal_id)
        #measure_GC_properties(gal_id)
    
    now4 = datetime.now()
    dt = (now4-now3)
    dt_string = now4.strftime("%d/%m/%Y %H:%M:%S")
    log.write('- STEP 4 DONE: '+dt_string+', TIME SPENT: '+str(dt.seconds)+'\n')

    finalize(gal_id)
    # step FINALE

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    log.write('END: '+dt_string+'\n\n')

    log.close()

if (MERGE_CATS == True):
    merge_cats()

if (MERGE_SIM_GC_CATS == True):
    merge_sims()

if (MERGE_GC_CATS == True):
    merge_gc_cats()
    #select_GC_candidadates_merged()
    #measure_GC_properties_merged()
