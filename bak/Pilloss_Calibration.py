import os, sys, glob
sys.path.insert( 1, os.getcwd() )

from __init__ import __version__
from utils import *
from operations import *

from input import params

def Calibration( work_path, slit_along, seeing, pixel_size ):
    '''
    '''
    
    # Wavelength Calibration
    # ======================
    wav_caled = wavelength_calibration( work_path, method = 'corr', inverse = True )
    
    # Extract Standard Spectrum
    # =========================
    std, _ =  load_image( os.path.join( work_path, 'corr/Std.1.fits' ), 0 )

    cnt_std_lst, cnt_std_err_lst = extract_oned_spec( imglist    = [std], 
                                                      seeing     = seeing / pixel_size, 
                                                      slit_along = slit_along, 
                                                      method     = 'auto', 
                                                      work_path  = work_path, 
                                                      imagetype  = 'Std'  )
    # Extract Object Spectrum
    # =======================
    objlist = list()
    for f in sorted( glob.glob( os.path.join( work_path, 'corr/Obj.*.fits' ) ) ):
        obj, _ =  load_image( f, 0 )
        objlist.append( obj )
    cnt_obj_lst, cnt_obj_err_lst = extract_oned_spec( imglist    = objlist, 
                                                      seeing     = seeing / pixel_size, 
                                                      slit_along = slit_along, 
                                                      method     = 'multi', 
                                                      work_path  = work_path, 
                                                      imagetype  = 'Obj'  )

    # Sensitivity Function
    # ====================
    sens = sens_func( wav_caled, cnt_std_lst[0], cnt_std_err_lst[0], inverse = True, exp = 600, airmass = 1.4 )

    # Flux Calibration
    # ================
    flux_calibration( wav_caled, cnt_obj_lst, cnt_obj_err_lst, sens, exp = 600, airmass = 1.2 )
    

def Pilloss_Calibration():

    Login_info( __version__ )

    Calibration( work_path  = os.getcwd(),
                 slit_along = params['slit_along'],
                 seeing     = 2.4,
                 pixel_size = 0.274,
               )



if __name__ == '__main__':

    Pilloss_Calibration()
