import os, argparse
from shutil import move

try:
    import numpy as np
except ImportError:
    print( 'Module `numpy` not found. Please install with: pip install numpy' )
    sys.exit()

try:
    from astropy.io import fits
except ImportError:
    print( 'Module `astropy` not found. Please install with: pip install astropy' )
    sys.exit()

from Pilloss_Prepare import check_file
    
def main( flist, idx, along, mode ):
    
    imglist = list()
    hdrlist = list()
    for f in flist:
        check_file( f )
        imglist.append( fits.getdata(   f, 0 ) )
        hdrlist.append( fits.getheader( f, 0 ) )
        move( f, 'bak/' )
        print( F'File `{f}` has been moved to `bak/`.' )
    
    if mode == 1:
        pass
    elif mode == 2:
        imglist = imglist[::-1]
        hdrlist = hdrlist[::-1]

    if along == 'row':
        data = np.hstack([ imglist[0][0, :, :idx], imglist[1][0, :, idx:] ])
        erro = np.hstack([ imglist[0][1, :, :idx], imglist[1][1, :, idx:] ])
    elif along == 'col':
        data = np.vstack([ imglist[0][0, :idx, :], imglist[1][0, idx:, :] ])
        erro = np.vstack([ imglist[0][1, :idx, :], imglist[1][1, idx:, :] ])

    fits.writeto( 'corr/Arc.1.fits', data = [data, erro], header = hdrlist[0], overwrite = True )
    print( F'The mosaic file `corr/Arc.1.fits` has been generated successfully.' )

if __name__ == '__main__':

    # Command line arguments
    parser = argparse.ArgumentParser( description = 'Concatenate two 1-D spectra at the pixel specified.' )
    parser.add_argument( '-index', help = 'Index of the point the two spectra files will be concatenated at.', default = '0' )
    parser.add_argument( '-along', help = '`row` or `col` along which the two spectra files will be concatenated', default = 'row' )
    parser.add_argument( '-mode', help = 'if `1`, FILE1[:index] will be concatenated with FILE2[index:]. If `2`, FILE2[:index] will be concatenated with FILE1[index:].', default = '1' )
    args = parser.parse_args()
    
    main( flist = ['corr/Arc.1.fits', 'corr/Arc.2.fits'], 
          idx   = int( args.index ), 
          along = args.along, 
          mode  = int( args.mode ) )