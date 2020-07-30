import os, sys

try:
    import numpy as np
except ImportError:
    print( 'Module `numpy` not found. Please install with: pip install numpy' )
    sys.exit()

try:
    from scipy.signal import find_peaks
    from scipy.optimize import curve_fit
    from scipy.ndimage import gaussian_filter
    from scipy.interpolate import interp1d
except ImportError:
    print( 'Module `scipy` not found. Please install with: pip install scipy' )
    sys.exit()

try:
    from prettytable import PrettyTable
    from prettytable import HEADER, NONE
except ImportError:
    print( 'Module `prettytable` not found. Please install with: pip install prettytable' )
    sys.exit()

try:
    from astropy.io import fits
    from astropy.time import Time
    from astropy.version import version
    from astropy.stats import sigma_clip
except ImportError:
    print( 'Module `astropy` not found. Please install with: pip install astropy' )
    sys.exit()

try:
    from spectres import spectres
except ImportError:
    print( 'Module `spectres` not found. Please install with: pip install spectres' )
    sys.exit()

try:
    import matplotlib.pyplot as plt
    # Set plot parameters
    plt.rcParams['axes.linewidth'] = 2.0
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'
except ImportError:
    print( 'Module `matplotlib` not found. Please install with: pip install matplotlib' )
    sys.exit()

try:
    from lacosmic import lacosmic
except ImportError:
    print( 'Module `lacosmic` not found. Please install with: pip install lacosmic' )
    sys.exit()

from utils import load_image, trim, imcmb, Gaussian, plot2d

def biascmb( work_path, blist, extension, trim_section ):

    '''
    A function to obtain combined bias.

    Parameters
    ----------
    work_path : string
        Path to the workspace.
    blist : array_like
        List of bias files.
    extension : integer or string
        The extension name in which data is stored.
    trim_section : list of tuples
        Specify column and row ranges.

    Returns
    -------
    bias : 3-D array
        Combined bias and its uncertainty.
    hdr : <class 'astropy.io.fits.header.Header'>
        Header.
    '''

    # Load bias
    # ---------
    bias = list( np.zeros([ 2, trim_section[1][1] - trim_section[1][0] + 1, trim_section[0][1] - trim_section[0][0] + 1 ]) for i in range( len(blist) ) )
    for i, b in enumerate(blist):
        print( F'\r[Bias Combination] Loading {i+1} of {len(blist)} bias files', end = '', flush = True )
        # Bias
        if i == 0:
            data,       hdr = load_image( os.path.join( work_path, F'{b}' ), extension )
            bias[i][0], hdr = trim( data, trim_section[0], trim_section[1], hdr )
        else:
            data,       _   = load_image( os.path.join( work_path, F'{b}' ), extension )
            bias[i][0], _   = trim( data, trim_section[0], trim_section[1] )
        # Uncertainty
        bias[i][1] = np.sqrt( bias[i][0] )
    # Print statistics
    tab = PrettyTable( hrules = HEADER, vrules = NONE )
    tab.field_names = [ 'FILE NAME', 'SHAPE', 'MEAN', 'STDDEV', 'MIN', 'MAX' ]
    for bname, b in zip( blist, bias ):
        tab.add_row([ bname, b[0].shape, round( b[0].mean(), 2 ), round( b[0].std(ddof=1), 2 ), b[0].min(), b[0].max() ])
    print( '\n\n' + tab.get_string() + '\n' )

    # Combine
    # -------
    print( '[Bias Combination] Combining' )
    bias, _, hdr = imcmb( bias, hdr, method = 'median', rejection = 'sigma_clip' )
    # Print statistics
    tab = PrettyTable( hrules = HEADER, vrules = NONE )
    tab.field_names = [ 'SHAPE', 'MEAN', 'STDDEV', 'MIN', 'MAX' ]
    tab.add_row([ bias[0].shape, round( bias[0].mean(), 2 ), round( bias[0].std(ddof=1), 2 ), round( bias[0].min(), 2 ), round( bias[0].max(), 2 ) ])
    print( '\n' + tab.get_string() + '\n' )

    # Writing to file
    # ---------------
    print( F'[Bias Combination] Writing to `{os.path.join( work_path, "bak/Bias.fits" )}`\n\n' )
    hdr['BITPIX'] = -32
    hdr['NAXIS']  = 3
    hdr['NAXIS3'] = 2
    hdr['EXTEND'] = ( False, 'File may contain extensions' )
    hdr['ORIGIN'] = ( F'Astropy v{version}', 'FITS file originator' )
    hdr['DATE'] = ( Time.now().fits, 'Date FITS file was generated' )
    hdr['COMMENT'] = 'Combined from '+ ', '.join(blist)
    del hdr['BSCALE'], hdr['BZERO']
    fits.writeto( os.path.join( work_path, 'bak/Bias.fits' ), data = bias, header = hdr, overwrite = True )

    return bias, hdr

def flatcmb( work_path, flist, extension, trim_section, bias ):

    '''
    A function to obtain combined flat-field.

    Parameters
    ----------
    work_path : string
        Path to the workspace.
    flist : array_like
        List of flat-field files.
    extension : integer or string
        The extension name in which data is stored.
    trim_section : list of tuples
        Specify column and row ranges.
    bias : 3-D array
        Combined bias which will be substracted from the flat-field
        and its uncertainty.

    Returns
    -------
    flat : 3-D array
        Combined flat-field and its uncertainty.
    '''

    # Load flat
    # ---------
    flat = list( np.zeros([ 2, bias.shape[1], bias.shape[2] ]) for i in range( len(flist) ) )
    for i, f in enumerate(flist):
        print( F'\r[Flat-field Combination] Loading {i+1} of {len(flist)} flat-field files', end = '', flush = True )
        if i == 0:
            data,       hdr = load_image( os.path.join( work_path, f ), extension )
            flat[i][0], hdr = trim( data, trim_section[0], trim_section[1], hdr )
        else:
            data,       _   = load_image( os.path.join( work_path, f ), extension )
            flat[i][0], _   = trim( data, trim_section[0], trim_section[1] )
        # Uncertainty
        flat[i][1] = np.sqrt( flat[i][0] )
    # Print statistics
    tab = PrettyTable( hrules = HEADER, vrules = NONE )
    tab.field_names = [ 'FILE NAME', 'SHAPE', 'MEAN', 'STDDEV', 'MIN', 'MAX' ]
    for fname, f in zip( flist, flat ):
        tab.add_row([ fname, f[0].shape, round( f[0].mean(), 2 ), round( f[0].std(ddof=1), 2 ), f[0].min(), f[0].max() ])
    print( '\n\n' + tab.get_string() + '\n' )

    # Combine
    # -------
    print( '[Flat-field Combination] Combining' )
    flat, _, hdr = imcmb( flat, hdr, method = 'mean', rejection = 'sigma_clip' )
    # Uncertainty
    flat[1] = np.sqrt( flat[1]**2 + bias[1]**2 )
    # Mean
    flat[0] = flat[0] - bias[0]
    # Print statistics
    tab = PrettyTable( hrules = HEADER, vrules = NONE )
    tab.field_names = [ 'SHAPE', 'MEAN', 'STDDEV', 'MIN', 'MAX' ]
    tab.add_row([ flat[0].shape, round( flat[0].mean(), 2 ), round( flat[0].std(ddof=1), 2 ), round( flat[0].min(), 2 ), round( flat[0].max(), 2 ) ])
    print( '\n' + tab.get_string() + '\n' )

    # Writing to file
    # ---------------
    print( F'[Flat-field Combination] Writing to `{os.path.join( work_path, "bak/Flat.fits" )}`\n\n' )
    hdr['BITPIX'] = -32
    hdr['NAXIS']  = 3
    hdr['NAXIS3'] = 2
    hdr['EXTEND'] = ( False, 'File may contain extensions' )
    hdr['ORIGIN'] = ( F'Astropy v{version}', 'FITS file originator' )
    hdr['DATE'] = ( Time.now().fits, 'Date FITS file was generated' )
    hdr['ZEROCOR'] = 'Zero level correction image is Bias.fits'
    hdr['COMMENT'] = 'Combined from '+ ', '.join(flist)
    del hdr['BSCALE'], hdr['BZERO']
    fits.writeto( os.path.join( work_path, 'bak/Flat.fits' ), data = flat, header = hdr, overwrite = True )

    return flat, hdr

def flatnorm( flat, sigma, work_path, hdr ):
    '''
    A function to normalize flatfield.

    Parameters
    ----------
    flat : 3-D array
        The flat to be normalized and its uncertainty.
    sigma : array_like
        sigma used in smoothing.    
    work_path : string
        Path to the workspace.
    hdr : <class 'astropy.io.fits.header.Header'>
        Header used to generate FITS file.

    Returns
    -------
    flat_norm : 3-D array
        Normalized flat and its uncertainty.
    '''

    print( '[Flat-field Normalization] Smoothing' )
    flat_smth = gaussian_filter( flat, sigma = ( 0, sigma[0], sigma[1] ), order = 0, mode = 'mirror' )
    flat_norm = flat[0] / flat_smth[0]
    # Uncertainty
    flat_norm_err = flat_norm * np.sqrt( ( flat[1] / flat[0] )**2 + ( flat_smth[1] / flat_smth[0] )**2 )
    flat_norm = np.stack([ flat_norm, flat_norm_err ])
    # Print statistics
    # ----------------
    tab = PrettyTable( hrules = HEADER, vrules = NONE )
    tab.field_names = [ 'SHAPE', 'MEAN', 'STDDEV', 'MIN', 'MAX' ]
    tab.add_row([ flat_norm[0].shape, round( flat_norm[0].mean(), 2 ), round( flat_norm[0].std(ddof=1), 2 ), round( flat_norm[0].min(), 2 ), round( flat_norm[0].max(), 2 ) ])
    print( '\n' + tab.get_string() + '\n' )

    print( F'[Flat-field Normalization] Writing to `{os.path.join( work_path, "bak/perFlat.fits" )}`\n\n' )
    hdr['BANDID1'] = F'Normalized flat - 2-D gaussian_filter, sigma = {sigma}'
    hdr['BANDID2'] = 'Uncertainty'
    fits.writeto( os.path.join( work_path, 'bak/perFlat.fits' ), data = flat_norm, header = hdr, overwrite = True )

    return flat_norm, flat_smth

def imcorr( work_path, flist, extension, imagetype, trim_section, bias, flat, slit_along, write = 0 ):
    '''
    A function to correct image.

    Parameters
    ----------
    work_path : string
        Path to the workspace.
    flist : array_like
        List of images.
    extension : integer or string
        The extension name in which data is stored.
    imagetype : string
        Type of the images.
    trim_section : list of tuples
        Specify column and row ranges.
    bias : 3-D array
        Combined bias which will be substracted from the flat-field
        and its uncertainty.
    flat : 3-D array
        Combined flat-field and its uncertainty.
    slit_along : `row` or `col`
        Specify direction of the slit.
    write : bool
        Flag to determine whether corrected images are written to files.

    Returns
    -------
    imglist : list
        List of corrected images and their uncertainty.
    hdrlist : list
        List of headers.
    '''

    # Load image
    # ----------
    imglist = list()
    hdrlist = list()
    for i, f in enumerate(flist):
        print( F'[Image Correction] {imagetype}: {os.path.join( work_path, f )}[{extension}]' )
        data, hdr = load_image( os.path.join( work_path, f ), extension )
        data, hdr = trim( data, trim_section[0], trim_section[1], hdr )
        # Corrected image
        img_corr = ( data - bias[0] ) / flat[0]
        # Uncertainty
        img_erro = np.sqrt( ( data + bias[1]**2 ) / ( data - bias[0] )**2 + ( flat[1] / flat[0] )**2 ) * np.abs( img_corr )
        imglist.append( np.stack([ img_corr, img_erro ]).astype( np.float32 ) )
        # Replace NaN with mean
        if np.isnan( imglist[i][1].sum() ):
            if slit_along == 'col': imglist[i][1] = np.where( np.isnan( imglist[i][1] ), np.ma.array( imglist[i][1], mask = np.isnan( imglist[i][1] ) ).mean( axis = 0 ), imglist[i][1] )
            if slit_along == 'row': imglist[i][1] = np.where( np.isnan( imglist[i][1] ), np.ma.array( imglist[i][1], mask = np.isnan( imglist[i][1] ) ).mean( axis = 1 )[:, np.newaxis], imglist[i][1] )
        hdrlist.append( hdr )
    # Print statistics
    tab = PrettyTable( hrules = HEADER, vrules = NONE )
    tab.field_names = [ 'FILE NAME', 'SHAPE', 'MEAN', 'STDDEV', 'MIN', 'MAX' ]
    for fname, img in zip( flist, imglist ):
        tab.add_row([ fname,         img[0].shape, round( img[0].mean(), 2 ), round( img[0].std(ddof=1), 2 ), img[0].min(), img[0].max() ])
        tab.add_row([ 'Uncertainty',         '\"', round( img[1].mean(), 2 ), round( img[1].std(ddof=1), 2 ), img[1].min(), img[1].max() ])
    print( '\n' + tab.get_string() + '\n' )

    # Edit header
    # -----------
    print( '[Image Correction] Edit header' )
    for hdr in hdrlist:
        hdr['BITPIX'] = -32
        hdr['NAXIS']  = 3
        hdr['NAXIS3'] = 2
        hdr['BANDID1'] = 'Preprocessed image - bias subtracted, flat-field corrected'
        hdr['BANDID2'] = 'Uncertainty'
        hdr['EXTEND'] = ( False, 'File may contain extensions' )
        hdr['ORIGIN'] = ( F'Astropy v{version}', 'FITS file originator' )
        hdr['DATE'] = ( Time.now().fits, 'Date FITS file was generated' )
        hdr['ZEROCOR'] = 'Zero level correction image is Bias.fits'
        hdr['FLATCOR'] = 'Flat-field correction image is perFlat.fits'
        del hdr['BSCALE'], hdr['BZERO']
    
    # Write to file
    # -------------
    if write:
        for i, f in enumerate(flist):
            fname = os.path.join( work_path, F'bak/{f[:-4]}_ppd.fits' )
            print( F'[Image Correction] Write to {fname}' )
            fits.writeto( fname, data = imglist[i], header = hdrlist[i], overwrite = True )
    
    print( '\n' )

    return imglist, hdrlist

def maskcr( imglist, hdrlist, imagetype, method, rdnoise, gain ):
    '''
    A function to mask cosmic ray pixels.

    Parameters
    ----------
    imglist : list
        List of images with cosmic ray pixels.
    hdrlist : list
        List of image headers.
    imagetype : string
        Type of the images.
    method : `Laplacian` or `median`
        If len(imglist) > 2, two methods are available, which are 
        Laplacian Edge Detection and median with sigma clipping
        denoted as `Laplacian` and `median` respectively. If 
        len(imglist) <= 2, Laplacian Edge Detection is used.
    rdnoise : float
        Read out noise.
    gain : float
        Photon gain.
    
    Returns
    -------
    data : 2-D array
        Mean object image.
    '''

    if len( imglist ) <= 2: method == 'Laplacian'
    
    if method == 'Laplacian':
        for i in range( len(imglist) ):
            print( F'[Remove Cosmic Rays] Process {i+1} of {len(imglist)} {imagetype} images using {method} method\n' )
            imglist[i][0], _ = lacosmic( imglist[i][0], 
                                         contrast           = 1, 
                                         cr_threshold       = 4.5, 
                                         neighbor_threshold = 0.5, 
                                         readnoise          = rdnoise, 
                                         effective_gain     = gain )
            hdrlist[i]['COMMENT'] = 'Cosmic ray pixels removed using Laplacian Edge Detection'
            print( '' )

    elif method == 'median':
        print( F'[Remove Cosmic Rays] Process {len(imglist)} {imagetype} images using {method} method\n' )
        data = np.array([ img[0] for img in imglist ])
        mask = np.zeros( data.shape, dtype = bool )
        for i in range( 3 ):
            # Median & standard deviation after maximum rejection
            maximum = data.max( axis = 0 )
            mask = mask | ( data == maximum )
            median = np.ma.median( np.ma.array( data, mask = mask ), axis = 0 )
            std    = np.ma.std(    np.ma.array( data, mask = mask ), axis = 0, ddof = 1 )
            # Sigma clipping & replace cosmic ray pixels with median value
            mask = ( data - median ) > 3 * std
            data_median = np.ma.median( np.ma.array( data, mask = mask ), axis = 0 )
            data[mask] = np.tile( data_median, [ data.shape[0], 1, 1 ] )[mask]
        for i in range( len(imglist) ):
            print( F'{i+1} of {len(imglist)} images: {mask[i].sum()} pixels masked [median with sigma clipping]\n' )
            imglist[i][0] = data[i]
            hdrlist[i]['COMMENT'] = 'Cosmic ray pixels removed using median with sigma clipping'
    
    else:
        raise ValueError( F'Unknown method `{method}`. Only `Laplacian` and `median` are available.' )

    print( '' )

    return imglist, hdrlist

def get_shift( imglist, slit_along, order, work_path ):
    '''
    A function to get zeropoint shift along dispersion axis.

    Parameters
    ----------
    imglist : list
        List of images used to determine shift curve.
    slit_along : `row` or `col`
        Specify direction of the slit.
    order : int
        Order used in polyfit.
    work_path : string
        Path to the workspace.
    
    Returns
    -------
    shift_fit : array_like
        Fitted zeropoint shift.
    '''
    
    shift = np.array([])
    for k, img in enumerate( imglist ):

        if slit_along == 'col':
            data = img[0]
        else:
            data = img[0].T
        # Find peaks
        # ----------
        print( F'[Zeropoint Shift] Find peaks on {k+1} of {len(imglist)} images' )
        peaks, properties = find_peaks( data.mean( axis = 0 ) / data.mean( axis = 0 ).max(), 
                                        height = ( 0.3, 0.8 ), 
                                        distance = data.shape[1] // 10,
                                        width = 0 )
        # Print peak info.
        tab = PrettyTable( hrules = HEADER, vrules = NONE )
        tab.field_names = [ 'CENTER', 'HEIGHT', 'WIDTH' ]
        for i, peak in enumerate( peaks ):
            tab.add_row([ peak, int( properties['peak_heights'][i] ), int( properties['widths'][i] ) ])
        print( '\n' + tab.get_string() + '\n' )
        
        # Derive zeropoint shift
        # ----------------------
        x = np.arange( data.shape[1] )
        mu = np.zeros([ data.shape[0], len( peaks ) ])
        for i in range( data.shape[0] ):
            # Peak center fitting
            print( F'\r[Zeropoint Shift] Peak center fitting on the ({i+1}/{data.shape[0]})-th row of {k+1}/{len(imglist)} images', end = '', flush = True )
            for j, peak in enumerate( peaks ):
                idxmin, idxmax = peak - int( properties['widths'][j] ), peak + int( properties['widths'][j] )
                popt, pcov = curve_fit( Gaussian, 
                                        x[idxmin:idxmax], 
                                        data[i, idxmin:idxmax], 
                                        bounds = ( [data[i, idxmin:idxmax].max()*0.5, x[idxmin:idxmax].min(), 0 ], 
                                                   [data[i, idxmin:idxmax].max()*1.5, x[idxmin:idxmax].max(), x[idxmin:idxmax].shape[0] ] ) )
                mu[i, j] = popt[1]

#                 if i == 260:
#                     plt.figure()
#                     plt.plot( x[idxmin:idxmax], data[i, idxmin:idxmax] )
#                     plt.plot( x[idxmin:idxmax], Gaussian(x[idxmin:idxmax], *popt), 'g--' )
#                     plt.savefig( './tmp1.png' )

        # Convert to relative shift
        shift = np.hstack([ shift.reshape( mu.shape[0], -1 ), ( mu - mu[mu.shape[0]//2] ) ])
        print( '' )
    # Mean
    shift_mean = shift.mean( axis = 1 )
    
    # Polyfit to shift curve
    # ----------------------
    print( F'[Zeropoint Shift] Fitting zeropoint shift iteratively (order = {order}, 5 iterations, nsigma = 1)' )
    y = np.arange( shift_mean.shape[0] ) + 1
    mask = np.ones( shift_mean.shape[0], dtype = bool )
    for k in range( 5 ):
        p = np.poly1d( np.polyfit( y[mask], shift_mean[mask], order ) )
        shift_fit = p( y )
        mask = mask & ( ~sigma_clip( shift_mean - shift_fit, sigma = 1, maxiters = 1, masked = True ).mask )
    
    np.savetxt( os.path.join( work_path, 'bak/zeropoint_shift.dat' ), shift_fit, fmt = '%15.8e' )

    # Plot
    # ----
    fig, ax = plt.subplots( 2, 1 )
    fig.subplots_adjust( hspace = 0 )

    for i in range( shift.shape[1] ):
        ax[0].plot( y, shift[:, i], 'r+' )
    ax[0].plot( y[~mask], shift_mean[~mask], '+', c = 'grey' )
    ax[0].plot( y[mask], shift_mean[mask], 'k+' )
    ax[0].plot( shift_fit, '-', c = 'yellow', lw = 2 )
    ax[0].set_xlim( y.min(), y.max() )
    ax[0].tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 14 )
    ax[0].set_xticklabels([])
    ax[0].set_ylabel( 'Displacement [px]', fontsize = 16 )

    ax[1].plot( y[mask],  shift_mean[mask]  - shift_fit[mask], 'k+' )
    ax[1].plot( y[~mask], shift_mean[~mask] - shift_fit[~mask], '+', c = 'grey' )
    ax[1].axhline( y = 0, ls = '--', c = 'yellow', lw = 2 )
    ax[1].set_xlim( y.min(), y.max() )
    ax[1].tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 14 )
    ax[1].set_xlabel( 'Slit', fontsize = 16 )
    ax[1].set_ylabel( 'Residuals [px]', fontsize = 16 )
    fig.align_ylabels()

    fig_path = os.path.join( work_path, 'figs' )
    print( F'[Zeropoint Shift] Plot zeropoint shift fitting to `{os.path.join( fig_path, "Zeropoint_shift_fitting.png" )}`\n\n' )
    if not os.path.exists( fig_path ):
        os.system( F'mkdir {fig_path}' )
    plt.savefig( os.path.join( fig_path, 'Zeropoint_shift_fitting.png' ), dpi = 144 )

    return shift_fit

def zp_shift_corr( imglist, hdrlist, imagetype, slit_along, shift, work_path ):
    '''
    '''

    for i, ( img, hdr ) in enumerate( zip( imglist, hdrlist ) ):

        print( F'[Zeropoint Shift] Correction for {i+1} of {len(imglist)} {imagetype} images' )
        if slit_along == 'row': 
            data = img[0].T; erro = img[1].T
        else:
            data = img[0]; erro = img[1]

        # Resample
        # --------
        new_data = np.zeros( data.shape )
        new_erro = np.zeros( erro.shape )
        x = np.arange( data.shape[1] )
        for j in range( data.shape[0] ):
            print( F'\r[Zeropoint Shift] Resampling ({j+1}/{data.shape[0]})-th row', end = '', flush = True )
            if shift[j] <= 0: fill = data[j][ 0]
            if shift[j] >  0: fill = data[j][-1]
            new_data[j], new_erro[j] = spectres( x, x - shift[j], data[j], erro[j], fill, verbose = False )
        print( '' )

        if slit_along == 'row': new_data = new_data.T; new_erro = new_erro.T
        new_data = ( new_data.astype( np.float32 ), new_erro.astype( np.float32 ) )

        # Write to file
        # -------------
        print( F'[Zeropoint Shift] Write to `{os.path.join( work_path, F"corr/{imagetype}.{str(i+1).zfill(len(str(len(imglist))))}.fits" )}`\n' )
        hdr['COMMENT'] = 'Zeropoint shifted'
        fits.writeto( os.path.join( work_path, F'corr/{imagetype}.{str(i+1).zfill(len(str(len(imglist))))}.fits' ), 
                      data      = new_data, 
                      header    = hdr, 
                      overwrite = True )

def wavelength_calibration( work_path, method = 'corr', inverse = True ):
    '''
    '''

    arc, _ = load_image( os.path.join( work_path, 'corr/Arc.1.fits' ), 0 )
    arc = np.median( arc[0], axis = 0 )
    arc = arc / arc.max()
    if inverse: arc = arc[::-1]
    
    if method == 'corr':

        wav, flx = np.loadtxt( './wavecal.dat' ).T

        x11 = np.arange( flx.shape[0] ) + 1
        x12 = np.arange( 1, flx.shape[0]+0.1, 0.1 )
        x21 = np.arange( arc.shape[0] ) + 1
        x22 = np.arange( 1, arc.shape[0]+0.1, 0.1 )
        flx_rsp = spectres( x12, x11, flx, verbose = False )
        arc_rsp = spectres( x22, x21, arc, verbose = False )

        idx = np.correlate( flx_rsp, arc_rsp, 'full' ).argmax() + 1
        wav_dense = interp1d( x11, wav, 'linear', bounds_error = False, fill_value = 'extrapolate' )( x12 )[np.max([0, idx-x22.shape[0]]):idx]
        if idx <= x12.shape[0]:
            idxmin = -wav_dense.shape[0]
            idxmax = x22.shape[0]
        elif idx > x12.shape[0]:
            if x11.shape[0] >= x21.shape[0]:
                idxmin = 0
                idxmax = wav_dense.shape[0]
            elif x11.shape[0] < x21.shape[0]:
                idxmin = -( idx - x12.shape[0] ) - wav_dense.shape[0]
                idxmax = -( idx - x12.shape[0] )
        wav_caled = interp1d( x22[idxmin:idxmax], wav_dense, 'linear', bounds_error = False, fill_value = 'extrapolate' )( x21 )
    
    fig, ax = plt.subplots( 1, 1, figsize = ( 12, 6 ) )
    ax.plot( wav_caled, arc )
    ax.plot( wav, flx )
#     ax.set_xticks( np.arange( 7500, 8550, 50 ), minor = True )
    ax.set_xlim( wav.min(), wav.max() )
    ax.tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 14 )
    ax.tick_params( which = 'minor', direction = 'in', top = True, right = True, length = 3, width = 1.5, labelsize = 14 )
    plt.savefig( './figs/tmp1.png' )

    if inverse: wav_caled = wav_caled[::-1]

    np.savetxt( os.path.join( work_path, 'bak/Wavelengths.dat' ), wav_caled, fmt = '%15.8e' )
    
    return wav_caled

def extract_oned_spec( imglist, seeing, slit_along, method, work_path, imagetype ):
    '''
    '''
    
    cnt_lst = list(); cnt_err_lst = list()
    for j, img in enumerate( imglist ):

        if slit_along == 'row': 
            data = img[0].T; erro = img[1].T
        else:
            data = img[0]; erro = img[1]

        if method == 'auto':

            print( F'[Extract 1-D Spectrum] Estimate background for {j+1} of {len(imglist)} {imagetype} image [{method}]' )
            guessidx = np.median( data, axis = 1 ).argmax()
            maxidx = np.zeros( data.shape[1] )
            for i in range( data.shape[1] ):
                idxmin, idxmax = ( guessidx - 2.5 * seeing ).astype( int ), ( guessidx + 2.5 * seeing ).astype( int )
                maxidx[i] = idxmin + data[idxmin:idxmax, i].argmax()

            idxbgmin = ( maxidx - 2.5 * seeing ).astype( int )
            idxbgmax = ( maxidx + 2.5 * seeing ).astype( int )

            x = np.arange( data.shape[0] ) + 1
            y = np.arange( data.shape[1] ) + 1
            cnt_tt     = np.zeros( data.shape[1] )
            cnt_tt_err = np.zeros( data.shape[1] )
            cnt_bg     = np.zeros( data.shape[1] )
            cnt_bg_err = np.zeros( data.shape[1] )
            data_bgsb  = np.zeros( data.shape    )
            for i in range( data.shape[1] ):
                idx = ( x < idxbgmin[i] ) | ( idxbgmax[i] < x )
                # Background fitting
                p = np.poly1d( np.polyfit( x[idx], data[idx, i], 1 ) )
                bg = p(x)
                # Background subtraction
                data_bgsb[:, i] = data[:, i] - bg
                # Background estimation
                cnt_bg[i] = bg[~idx].sum()
                # Background uncertainty
                cnt_bg_err[i] = np.median( erro[idx, i] ) * np.sqrt( ( ~idx ).sum() )
                # Total flux estimation
                cnt_tt[i] = data[~idx, i].sum()
                # Total flux uncertainty
                cnt_tt_err[i] = np.sqrt( ( erro[~idx, i]**2 ).sum() )
            
            print( F'[Extract 1-D Spectrum] Estimate target flux and uncertainty for {j+1} of {len(imglist)} {imagetype} image [{method}]' )
            # Source flux estimation
            cnt_sc = cnt_tt - cnt_bg
            # Source flux uncertainty
            cnt_sc_err = np.sqrt( cnt_tt_err**2 + cnt_bg_err**2 )

        elif method == 'multi':

            print( F'[Extract 1-D Spectrum] Estimate background for {j+1} of {len(imglist)} {imagetype} image [{method}]' )
            guessidx = np.median( data, axis = 1 ).argmax()
            maxidx = np.zeros( data.shape[1] )
            for i in range( data.shape[1] ):
                idxmin, idxmax = ( guessidx - 2.5 * seeing ).astype( int ), ( guessidx + 2.5 * seeing ).astype( int )
                maxidx[i] = idxmin + data[idxmin:idxmax, i].argmax()

            idxbgmin = ( maxidx - 60 * seeing ).astype( int )
            idxbgmax = ( maxidx + 60 * seeing ).astype( int )

            idxscmin = ( maxidx - 10 * seeing ).astype( int )
            idxscmax = ( maxidx + 10 * seeing ).astype( int )

            x = np.arange( data.shape[0] ) + 1
            y = np.arange( data.shape[1] ) + 1
            cnt_tt     = np.zeros( data.shape[1] )
            cnt_tt_err = np.zeros( data.shape[1] )
            cnt_bg     = np.zeros( data.shape[1] )
            cnt_bg_err = np.zeros( data.shape[1] )
            data_bgsb  = np.zeros( data.shape    )
            for i in range( data.shape[1] ):
                idxbg = ( x < idxbgmin[i] ) | ( idxbgmax[i] < x )
                idxsc = ( x < idxscmin[i] ) | ( idxscmax[i] < x )
                # Background fitting
                p = np.poly1d( np.polyfit( x[idxbg], data[idxbg, i], 1 ) )
                bg = p(x)
                # Background subtraction
                data_bgsb[:, i] = data[:, i] - bg
                # Background estimation
                cnt_bg[i] = bg[~idxsc].sum()
                # Background uncertainty
                cnt_bg_err[i] = np.median( erro[idxbg, i] ) * np.sqrt( ( ~idxsc ).sum() )
                # Total flux estimation
                cnt_tt[i] = data[~idxsc, i].sum()
                # Total flux uncertainty
                cnt_tt_err[i] = np.sqrt( ( erro[~idxsc, i]**2 ).sum() )
            
            print( F'[Extract 1-D Spectrum] Estimate target flux and uncertainty for {j+1} of {len(imglist)} {imagetype} image [{method}]' )
            # Source flux estimation
            cnt_sc = cnt_tt - cnt_bg
            # Source flux uncertainty
            cnt_sc_err = np.sqrt( cnt_tt_err**2 + cnt_bg_err**2 )
        else:
            raise ValueError( 'Unknown method. Only `auto` and `multi` are available.' )
        
        # Write to file
        file_path = os.path.join( work_path, F'bak/1-D_Extracted_{imagetype}_Spectrum_{str(j+1).zfill(len(str(len(imglist))))}.dat' )
        print( F'[Extract 1-D Spectrum] Write 1-D spectrum extracted from {j+1} of {len(imglist)} {imagetype} image to `{file_path}`' )
        np.savetxt( file_path, np.vstack([ cnt_sc, cnt_sc_err, cnt_bg, cnt_bg_err ]).T, fmt = '%15.8e' )
        
        # Plots
        fig_path  = os.path.join( work_path, 'figs' )
        # 1-D plot
        title = F'1-D Extracted {imagetype} Spectrum {str(j+1).zfill(len(str(len(imglist))))}'
        print( F'[Extract 1-D Spectrum] Plot 1-D spectrum extracted from {j+1} of {len(imglist)} {imagetype} image to `{os.path.join( fig_path, F"{title}.png".replace( " ", "_" ) )}`\n\n' )
        fig, ax = plt.subplots( 1, 1, figsize = ( 8, 5 ) )
        ax.step( y, cnt_tt, 'k', where = 'mid', label = 'Total' )
        ax.step( y, cnt_sc, where = 'mid', label = 'Target' )
        ax.fill_between( y, cnt_sc - cnt_sc_err, cnt_sc + cnt_sc_err, alpha = 0.2 )
        ax.step( y, cnt_bg, where = 'mid', label = 'Background' )
        ax.grid( axis = 'both', color = '0.95' )
        ax.set_xlim( y.min(), y.max() )
        ax.tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 14 )
        ax.set_xlabel( 'Dispersion Direction [px]', fontsize = 16 ); ax.set_ylabel( 'Counts', fontsize = 16 )
        ax.legend( fontsize = 14 )
        ax.set_title( F'Extracted {imagetype} Spectrum', fontsize = 18 )
        plt.savefig( os.path.join( fig_path, F'{title}.png'.replace( ' ', '_' ) ), dpi = 144 )
        plt.close()
        
        # 2-D plot
        title = F'2-D Background Subtracted {imagetype} image {str(j+1).zfill(len(str(len(imglist))))}'
        print( F'[Extract 1-D Spectrum] Plot 2-D spectrum extracted from {j+1} of {len(imglist)} {imagetype} image to `{os.path.join( fig_path, F"{title}.png".replace( " ", "_" ) )}`\n\n' )
        fig, ax = plt.subplots( 1, 1 )
        fig.subplots_adjust( right = 0.8 )

        # Image
        im = ax.imshow( data_bgsb, cmap = 'Greys_r', origin = 'lower', extent = ( 0.5, data_bgsb.shape[1]+0.5, 0.5, data_bgsb.shape[0]+0.5) )
        # Background aperture
        ax.plot( np.arange( data_bgsb.shape[1] ) + 1, idxbgmin, 'lightskyblue', ls = '--' )
        ax.plot( np.arange( data_bgsb.shape[1] ) + 1, idxbgmax, 'lightskyblue', ls = '--' )
        # Target aperture
        if method == 'auto': 
            ax.plot( np.arange( data_bgsb.shape[1] ) + 1, idxbgmin, 'yellow', ls = '-' )
            ax.plot( np.arange( data_bgsb.shape[1] ) + 1, idxbgmax, 'yellow', ls = '-' )
        else:
            ax.plot( np.arange( data_bgsb.shape[1] ) + 1, idxscmin, 'yellow', ls = '-' )
            ax.plot( np.arange( data_bgsb.shape[1] ) + 1, idxscmax, 'yellow', ls = '-' )
        # Colorbar
        cax = fig.add_axes([ ax.get_position().x1 + 0.02, ax.get_position().y0, 0.04, ax.get_position().height ])
        cb = fig.colorbar( im, cax = cax )
        # Settings
        ax.tick_params( which = 'major', direction = 'in', color = 'w', top = True, right = True, length = 7, width = 1.5, labelsize = 14 )
        if slit_along == 'col': ax.set_xlabel( 'Dispersion', fontsize = 18 ); ax.set_ylabel( 'Slit', fontsize = 18 )
        if slit_along == 'row': ax.set_xlabel( 'Slit', fontsize = 18 ); ax.set_ylabel( 'Dispersion', fontsize = 18 )
        ax.set_title( title, fontsize = 18 )
        cb.ax.tick_params( which = 'major', direction = 'in', color = 'w', right = True, length = 7, width = 1.5, labelsize = 14 )
        cb.ax.set_ylabel( 'ADU', fontsize = 18 )
        plt.savefig( os.path.join( fig_path, F'{title}.png'.replace( ' ', '_' ) ), dpi = 144 )

        cnt_lst.append( cnt_sc )
        cnt_err_lst.append( cnt_sc_err )

    return cnt_lst, cnt_err_lst

def sens_func( wav, cnt, cnt_err, inverse, exp, airmass = 1.4 ):
    '''
    '''

    # Speed of light
    c = 2.99792458e+10 # [cm/s]
    
    if inverse: wav = wav[::-1]; cnt = cnt[::-1]; cnt_err = cnt_err[::-1]

    # Atmospheric extinction correction
    # ---------------------------------
    wav_ext, ext = np.loadtxt( './baoextinct.dat' ).T
    ext = interp1d( wav_ext, ext, kind = 'quadratic' )( wav )
    ext_corr_factor = 10**( 0.4 * airmass * ext )
    cnt = cnt * ext_corr_factor
    
    # Bandpass
    # --------
    dwav = np.abs( np.diff( wav ) )
    dwav = np.hstack([ dwav[0], dwav, dwav[-1] ])
    dwav = ( dwav[:-1] + dwav[1:] ) / 2
    
    cnt = cnt / dwav / exp # [counts/A/s]

    # Sensitivity function
    # --------------------
    wav_mag, mag, bp = np.loadtxt( './he3.dat', skiprows = 1 ).T
    flx_mod = 10**( -0.4 * mag ) * 3631e-23 * c / wav_mag**2  * 1e8 # [erg/cm2/s/A]
    flx_obs = np.zeros( flx_mod.shape[0] )
    bins = np.vstack([ wav_mag - bp//2, wav_mag + bp//2 ]).T
    for i, bin in enumerate( bins ):
        idx = ( bin[0] < wav ) & ( wav < bin[1] )
        edges = interp1d( wav, cnt, 'linear', bounds_error = False, fill_value = np.nan )( bin )
        flx_obs[i] = np.trapz( np.hstack([ edges[0], cnt[idx], edges[1] ]), x = np.hstack([ bin[0], wav[idx], bin[1] ]) ) / ( bin[1] - bin[0] )
    sen = flx_obs / flx_mod
    idx = np.isnan( sen )
    wav_sen = wav_mag[~idx]; sen = 2.5 * np.log10( sen[~idx] )
    
    mask = np.ones( wav_sen.shape[0], dtype = bool )
    for i in range( 5 ):
        p = np.poly1d( np.polyfit( wav_sen[mask], sen[mask], 12 ) )
        mask = mask & ~sigma_clip( sen - p( wav_sen ), sigma = 2, maxiters = 1, masked = True ).mask
    sen_fit = p( wav )

    fig, ax = plt.subplots( 2, 1, figsize = (8, 8) )

    ax[0].plot( wav_sen[ mask], sen[ mask], 'b+', ms = 10 )
    ax[0].plot( wav_sen[~mask], sen[~mask], 'r+', ms = 10 )
    ax[0].plot( wav, sen_fit, 'k-', zorder = 0 )
    ax[0].grid( axis = 'both', color = '0.95' )
    ax[0].set_xlim( wav.min(), wav.max() )
    ax[0].tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 14 )
    ax[0].set_ylabel( 'Sensitivity Function', fontsize = 16 )

    ax[1].plot( wav_sen[ mask], ( sen - p( wav_sen ) )[ mask], 'b+', ms = 10 )
    ax[1].plot( wav_sen[~mask], ( sen - p( wav_sen ) )[~mask], 'r+', ms = 10 )
    ax[1].axhline( y = 0, c = 'grey', ls = '--' , zorder = 0 )
    ax[1].grid( axis = 'both', color = '0.95' )
    ax[1].set_xlim( wav.min(), wav.max() )
    ax[1].set_ylim( ( sen - p( wav_sen ) )[mask].min() * 0.8, ( sen - p( wav_sen ) )[mask].max() * 1.2 )
    ax[1].tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 14 )
    ax[1].set_xlabel( 'Wavelength [$\\mathrm{\\AA}$]', fontsize = 16 ); 
    ax[1].set_ylabel( 'Residuals', fontsize = 16 )
    fig.align_ylabels()

    plt.savefig( './figs/tmp.png' )
    
    if inverse: sen_fit = sen_fit[::-1]
    
    return sen_fit

def flux_calibration( wav, cnt_lst, cnt_err_lst, sens, exp, airmass = 1.4 ):
    '''
    '''
    
    for i, ( cnt, cnt_err ) in enumerate( zip( cnt_lst, cnt_err_lst ) ):

        # Atmospheric extinction correction
        # ---------------------------------
        wav_ext, ext = np.loadtxt( './baoextinct.dat' ).T
        ext = interp1d( wav_ext, ext, kind = 'quadratic' )( wav )
        ext_corr_factor = 10**( 0.4 * airmass * ext )
        cnt = cnt * ext_corr_factor

        # Bandpass
        # --------
        dwav = np.abs( np.diff( wav ) )
        dwav = np.hstack([ dwav[0], dwav, dwav[-1] ])
        dwav = ( dwav[:-1] + dwav[1:] ) / 2

        # Calibration
        # -----------
        cnt = cnt / dwav / exp # [counts/A/s]
        flx = cnt * 10**( -0.4 * sens )

        idx = ( 3700 < wav ) & ( wav < 9000 )
        fig, ax = plt.subplots( 1, 1, figsize = (8, 4) )
        ax.step( wav[idx], flx[idx], 'k-', where = 'mid' )
        ax.set_ylim( 0, 1e-14 )
        plt.savefig( F'./figs/tmp2.{i+1}.png' )

        np.savetxt( F'./tmp.{i+1}.dat', np.vstack([wav[idx], flx[idx]]).T )