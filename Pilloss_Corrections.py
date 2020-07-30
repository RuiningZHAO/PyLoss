import os, sys, argparse
sys.path.insert( 1, os.getcwd() )

from __init__ import __version__
from utils import *
from operations import *

from input import params

def Corrections( work_path, lists, ext, trim_section, slit_along, rdnoise, gain, correction = True ):
    '''
    '''

    # BIAS
    # ====
    # 1. Combine
    bias, bias_hdr = biascmb( work_path, lists['bias'], ext, trim_section )
    # 2. Plot
    plot2d( bias[0], work_path, 'Bias', slit_along, show = 0, save = 1 )

    # FLAT-FIELD
    # ==========
    # 1. Combine
    flat, flat_hdr = flatcmb( work_path, lists['flat'], ext, trim_section, bias )
    plot2d( flat[0], work_path, 'Flat-field', slit_along, show = 0, save = 1 )
    # 2. Normalize
    flat_norm, flat_smth = flatnorm( flat, (20, 20), work_path, flat_hdr )
    # 3. Plot
    plot2d( flat_smth[0], work_path, 'Smoothed Flat-field',   slit_along, show = 0, save = 1 )
    plot2d( flat_norm[0], work_path, 'Normalized Flat-field', slit_along, show = 0, save = 1 )

    # CORRECTIONS
    # ===========
    if correction:
        # 1. Load preprocessed spectra
        arcs, arc_hdrs = imcorr( work_path, lists['arc'], ext, 'Arc', trim_section, bias, flat_norm, slit_along, write = 0 )
        objs, obj_hdrs = imcorr( work_path, lists['obj'], ext, 'Obj', trim_section, bias, flat_norm, slit_along, write = 0 )
        stds, std_hdrs = imcorr( work_path, lists['std'], ext, 'Std', trim_section, bias, flat_norm, slit_along, write = 0 )
        # 2. Mask cosmic rays
        objs, obj_hdrs = maskcr( objs, obj_hdrs, 'Obj', method = 'Laplacian', rdnoise = rdnoise, gain = gain )
        stds, std_hdrs = maskcr( stds, std_hdrs, 'Std', method = 'Laplacian', rdnoise = rdnoise, gain = gain )
        # 3. Get fitted zeropoint shift
        shift = get_shift( arcs, slit_along, 7, work_path )
        # 4. Zeropoint shift correction
        zp_shift_corr( arcs, arc_hdrs, 'Arc', slit_along, shift, work_path )
        zp_shift_corr( objs, obj_hdrs, 'Obj', slit_along, shift, work_path )
        zp_shift_corr( stds, std_hdrs, 'Std', slit_along, shift, work_path )

def Pilloss_Corrections( correction = True ):

    Login_info( __version__ )
    
    lists = load_lists( work_path = os.getcwd(), listnames = ['bias', 'flat', 'arc', 'obj', 'std'] )

    Corrections( work_path    = os.getcwd(), 
                 lists        = lists, 
                 ext          = params['extension'], 
                 trim_section = params['trim_section'], 
                 slit_along   = params['slit_along'],
                 rdnoise      = params['rdnoise'],
                 gain         = params['gain'],
                 correction   = correction )

if __name__ == '__main__':

    # Command line arguments
    parser = argparse.ArgumentParser( description = 'Correction in Long-slit spectroscopy' )
    parser.add_argument( '-correction', help = 'Do correction. If `0`, only correction files generated.', default = 1 )

    args = parser.parse_args()
    correction = int( args.correction )

    Pilloss_Corrections( correction )
