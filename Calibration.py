import os, sys
from glob import glob
from copy import deepcopy

try:
    import numpy as np
except ImportError:
    print( 'Module `numpy` not found. Please install with: pip install numpy' )
    sys.exit()

try:
    from scipy.signal import find_peaks
    from scipy.optimize import curve_fit
    from scipy.ndimage import gaussian_filter
    from scipy.interpolate import interp1d, make_lsq_spline
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
    from astropy.stats import sigma_clip
    from astropy.modeling.models import Gaussian2D
    from astropy.stats import gaussian_fwhm_to_sigma
except ImportError:
    print( 'Module `astropy` not found. Please install with: pip install astropy' )
    sys.exit()
    
try:
    from spectres import spectres
except ImportError:
    print( 'Module `spectres` not found. Please install with: pip install spectres' )
    sys.exit()

try:
    from photutils import RectangularAperture, aperture_photometry
except ImportError:
    print( 'Module `photutils` not found. Please install with: pip install photutils' )
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

from .utils import load_image, imcmb, Gaussian, plot2d

def get_slit_loss( fwhm, sw ):
    '''
    '''

    sig = fwhm * gaussian_fwhm_to_sigma

    ny = nx = np.int( 10 * sig )
    Y, X = np.mgrid[-ny:ny+1, -nx-10:nx+1+10]

    g = Gaussian2D( amplitude = 1 / ( 2 * np.pi * sig**2 ), 
                    x_mean    = 0, 
                    y_mean    = 0, 
                    x_stddev  = sig, 
                    y_stddev  = sig,
                    theta     = 0 )
    data = g( X, Y )

    ( y, x ) = np.unravel_index( np.argmax( data, axis = None ), data.shape )

    aper = RectangularAperture( ( x, y ), w = sw, h = data.shape[0] )

    # HM = data[ind] / 2
    # FW = np.where( HM <= data[ind[0]] )[0][-1] - np.where( HM <= data[ind[0]] )[0][ 0]
    # print( FW * gaussian_fwhm_to_sigma / sig )

    # Aperture photometry
    phot_table = aperture_photometry( data, aper, method = 'subpixel', subpixels = 100 )

    slit_loss =  data.sum() - phot_table['aperture_sum'][0]

#     fig, ax = plt.subplots( 1, 1, figsize = ( 10, 10 ) )
#     ax.imshow( data )
#     aper.plot( axes = ax, color = 'red', lw = 2 )
#     plt.show()

    return slit_loss

def extract_oned_spec( path, imgtype, slit_along, method = 'aper', aper_param = None, seeing = None, order = None, show = True, save = True ):
    '''
    '''
    
    filelist = sorted( glob( F'{path}/{imgtype}*.fits' ) )
    
    if len( filelist ) == 0: raise FileExistsError( F'File {F"{path}/{imgtype}*.fits"} do not exist.' )

    cnts_sc = list(); errs_sc = list()
    cnts_bg = list(); errs_bg = list()
    for i, f in enumerate( filelist ):

        # Load image
        ( img, err ), _ = load_image( f, 0 )

        # Transpose
        if slit_along == 'row': img = img.T; err = err.T
        
        x = np.arange( img.shape[1] ); y = np.arange( img.shape[0] )
        
        if method == 'median':
            print( F'[Extract 1-D Spectrum] Extract {i+1} of {len(filelist)} {imgtype} image [{method}]' )
            cnt_sc = np.median( img, axis = 0 )
            err_sc = np.sqrt( ( err**2 ).sum( axis = 0 ) ) / err.shape[0]

        elif method == 'aper':
            print( F'[Extract 1-D Spectrum] Estimate background for {i+1} of {len(filelist)} {imgtype} image [{method}]' )
            naper, aper_sc, aper_bg = aper_param
            guessidx = np.median( img, axis = 1 ).argmax()
            maxidx = np.zeros( img.shape[1] )
            for j in range( img.shape[1] ):
                idxmin, idxmax = ( guessidx - 1.5 * seeing ).astype( int ), ( guessidx + 1.5 * seeing ).astype( int )
                maxidx[j] = idxmin + gaussian_filter( img[idxmin:idxmax, j], sigma = seeing ).argmax()

#             p = np.poly1d( np.polyfit( x, maxidx, 3 ) )
            knots = np.r_[ (x[0], ) * ( order + 1 ), ( x[-1], ) * ( order + 1 ) ]
            spl = make_lsq_spline( x, maxidx, t = knots, k = order )
            maxloc = spl(x)

            locbgmin = ( maxloc - aper_bg/2 )
            locbgmax = ( maxloc + aper_bg/2 )

            cnt_tt = np.zeros([ naper, img.shape[1] ]); err_tt = np.zeros([ naper, img.shape[1] ])
            cnt_bg = np.zeros([ naper, img.shape[1] ]); err_bg = np.zeros([ naper, img.shape[1] ])
            img_bgsb = np.zeros( img.shape )
            apers = np.zeros([ naper+1, img.shape[1] ])
            for j in range( img.shape[1] ):
                # Background fitting
                idxbg = ( y < locbgmin[j] ) | ( locbgmax[j] < y )
                p = np.poly1d( np.polyfit( y[idxbg], img[idxbg, j], 1 ) )
                bg = p(y)
                # Extract 1-D spectra
                apers[:, j] = np.linspace( -aper_sc/2, aper_sc/2, naper + 1 ) + maxloc[j]
                locscmins, locscmaxs = apers[:-1, j], apers[1:, j]
                for k, ( locscmin, locscmax ) in enumerate( zip( locscmins, locscmaxs ) ):
                    idxsc = ( locscmin < y ) & ( y < locscmax )
                    # Background estimation
                    cnt_bg[k, j] = bg[idxsc].sum() + \
                                   bg[y[idxsc][ 0]-1] * ( y[idxsc][ 0] - locscmin ) + \
                                   bg[y[idxsc][-1]+1] * ( locscmax - y[idxsc][-1] )
                    # Background uncertainty
                    err_bg[k, j] = np.median( err[idxbg, j] ) * np.sqrt( locscmax - locscmin + 1 )
                    # Total flux estimation
                    cnt_tt[k, j] = img[idxsc, j].sum() + \
                                   img[y[idxsc][ 0]-1, j] * ( y[idxsc][ 0] - locscmin ) +\
                                   img[y[idxsc][-1]+1, j] * ( locscmax - y[idxsc][-1] )
                    # Total flux uncertainty
                    err_tt[k, j] = np.sqrt( ( err[idxsc, j]**2 ).sum() + \
                                            ( err[y[idxsc][ 0]-1, j] * ( y[idxsc][ 0] - locscmin ) )**2 + \
                                            ( err[y[idxsc][-1]+1, j] * ( locscmax - y[idxsc][-1] ) )**2 )
                # Background subtraction
                img_bgsb[:, j] = img[:, j] - bg
            
            print( F'[Extract 1-D Spectrum] Estimate target flux and uncertainty for {i+1} of {len(filelist)} {imgtype} image [{method}]' )
            # Source flux estimation
            cnt_sc = cnt_tt - cnt_bg
            # Source flux uncertainty
            err_sc = np.sqrt( err_tt**2 + err_bg**2 )

            cnts_bg.append( cnt_bg ); errs_bg.append( err_bg )

            # 2-D plot
            title = F'2-D Background Subtracted {imgtype} image {str(i+1).zfill(len(str(len(filelist))))}'
            print( F'[Extract 1-D Spectrum] Plot 2-D background subtracted {imgtype} image from {i+1} of {len(filelist)} {imgtype} image', end = '' )

            fig, ax = plt.subplots( 1, 1, figsize = ( 10, 8 ) )
            fig.subplots_adjust( right = 0.8 )
            # Image
            im = ax.imshow( img_bgsb, cmap = 'Greys_r', origin = 'lower', extent = ( 0.5, img_bgsb.shape[1]+0.5, 0.5, img_bgsb.shape[0]+0.5) )
            # Background aperture
            ax.plot( x+1, locbgmin+1, 'lightskyblue', ls = '--' ); ax.plot( x+1, locbgmax+1, 'lightskyblue', ls = '--' )
            # Target aperture
            for aper in apers:
                ax.plot( x+1, aper+1, 'yellow', ls = '-' )
            # Colorbar
            cax = fig.add_axes([ ax.get_position().x1 + 0.02, ax.get_position().y0, 0.04, ax.get_position().height ])
            cb = fig.colorbar( im, cax = cax )
            # Settings
            ax.tick_params( which = 'major', direction = 'in', color = 'w', top = True, right = True, length = 7, width = 1.5, labelsize = 18 )
            if slit_along == 'col': ax.set_xlabel( 'Dispersion', fontsize = 22 ); ax.set_ylabel( 'Slit', fontsize = 22 )
            if slit_along == 'row': ax.set_xlabel( 'Slit', fontsize = 22 ); ax.set_ylabel( 'Dispersion', fontsize = 22 )
            ax.set_title( title, fontsize = 24 )
            cb.ax.tick_params( which = 'major', direction = 'in', color = 'w', right = True, length = 7, width = 1.5, labelsize = 18 )
            cb.ax.set_ylabel( 'ADU', fontsize = 22 )
            if save:
                print( F' to `{ F"figs/{title}.png".replace( " ", "_" ) }`' )
                if not os.path.exists( 'figs' ): os.makedirs( 'figs' )
                plt.savefig( F'figs/{title}.png'.replace( ' ', '_' ), dpi = 144 )
            else:
                print( '' )
            if show:
                print( F'[Extract 1-D Spectrum] Show plot' )
                plt.show()
            plt.close()

        else:
            raise ValueError( 'Unknown method. Only methods in [`median`, `aper`] are available.' )

        cnts_sc.append( cnt_sc ); errs_sc.append( err_sc )

        # 1-D plot
        title = F'1-D Extracted {imgtype} Spectrum {str(i+1).zfill(len(str(len(filelist))))}'
        print( F'[Extract 1-D Spectrum] Plot 1-D spectrum extracted from {i+1} of {len(filelist)} {imgtype} image', end = '' )

        fig, ax = plt.subplots( 1, 1, figsize = ( 12, 6 ) )
        if method == 'median':
            # Target
            ax.step( x+1, cnt_sc, where = 'mid', label = 'Target', zorder = 3 )
            ax.fill_between( x+1, cnt_sc - err_sc, cnt_sc + err_sc, color = 'C0', alpha = 0.2, zorder = 2 )
        else:
            for m in range( cnt_sc.shape[0] ):
                # Target
                if m == 0:
                    ax.step( x+1, cnt_sc[m], where = 'mid', color = 'C0', zorder = 3, label = 'Target' )
                else:
                    ax.step( x+1, cnt_sc[m], where = 'mid', color = 'C0', zorder = 3 )
                ax.fill_between( x+1, cnt_sc[m] - err_sc[m], cnt_sc[m] + err_sc[m], color = 'C0', alpha = 0.2, zorder = 2 )
            # Background
            ax.step( x+1, np.median( cnt_bg, axis = 0 ), where = 'mid', color = 'C1', label = 'Background', zorder = 1 )
        # Settings
        ax.grid( axis = 'both', color = '0.95', zorder = 0 )
        ax.set_xlim( x.min()+1, x.max()+1 )
        ax.tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 18 )
        ax.set_xlabel( 'Dispersion Direction [px]', fontsize = 22 ); ax.set_ylabel( 'Counts', fontsize = 22 )
        ax.legend( fontsize = 22 )
        ax.set_title( title, fontsize = 24 )
        
        if save:
            print( F' to `{ F"figs/{title}.png".replace( " ", "_" ) }`' )
            if not os.path.exists( 'figs' ): os.makedirs( 'figs' )
            plt.savefig( F'figs/{title}.png'.replace( ' ', '_' ), dpi = 144 )
        else:
            print( '' )
        if show:
            print( '[Extract 1-D Spectrum] Show plot' )
            plt.show()
        plt.close()

        # Write to file
        if not os.path.exists( 'bak' ): os.makedirs( 'bak' )
        file_path = F'bak/1-D_Extracted_{imgtype}_Spectrum_{str(i+1).zfill(len(str(len(filelist))))}.dat'
        print( F'[Extract 1-D Spectrum] Write 1-D spectrum extracted from {i+1} of {len(filelist)} {imgtype} image to `{file_path}`' )
        if method == 'median':
            np.savetxt( file_path, np.vstack([ cnt_sc, err_sc ]).T, fmt = '%15.8e' )
        else:
            np.savetxt( file_path, np.vstack([ cnt_sc, err_sc, cnt_bg, err_bg ]).T, fmt = '%15.8e' )
    
    if method == 'median':
        return ( cnts_sc, errs_sc )
    else:
        return ( cnts_sc, errs_sc ), ( cnts_bg, errs_bg )

def wavelength_calibration( arc, method, wavfile, precision = None, inverse = False, show = 1, save = 0 ):
    '''
    '''

    if inverse: arc = arc[::-1]

    # Normalize
    arc = arc / arc.max()
    # Calibrate
    if method == 'corr':
        print( '[Wavelength calibration] Using cross correlation algorithm' )
        print( '[Wavelength calibration] Load calibrated wavelength grid' )
        wav, flx = np.loadtxt( wavfile ).T
        print( F'[Wavelength calibration] Resample to denser grid [{precision} px]' )
        x11 = np.arange( flx.shape[0] ) + 1
        x12 = np.arange( 1, flx.shape[0] + precision, precision )  # dense
        x21 = np.arange( arc.shape[0] ) + 1
        x22 = np.arange( 1, arc.shape[0] + precision, precision )  # dense
        flx_rsp = spectres( x12, x11, flx, verbose = False )
        arc_rsp = spectres( x22, x21, arc, verbose = False )
        print( '[Wavelength calibration] Locate maximum' )
        idx = np.correlate( flx_rsp, arc_rsp, 'full' ).argmax() + 1
        print( '[Wavelength calibration] Interpolate back to the origin grid' )
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

        # Plot
        print( '[Wavelength calibration] Plot wavelength calibrated 1-D spectrum', end = '' )
        fig, ax = plt.subplots( 1, 1, figsize = ( 12, 6 ) )
        # Newly calibrated
        ax.plot( wav_caled, arc, color = 'C0', zorder = 3, label = 'Newly calibrated' )
        # Archived
        ax.plot( wav, flx, color = 'C1', zorder = 2, label = 'Archived' )
        # Settings
        ax.grid( axis = 'both', color = '0.95', zorder = 0 )
        ax.set_xlim( wav.min(), wav.max() )
        ax.tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 18 )
        ax.set_xlabel( 'Wavelength [$\\mathrm{\\AA}$]', fontsize = 22 )
        ax.set_ylabel( 'Relative counts', fontsize = 22 )
        ax.legend( fontsize = 22 )
        ax.set_title( 'Wavelength Calibrated 1-D Spectrum', fontsize = 24 )
        if save:
            print( ' to `figs/Wavelength_Calibrated_1-D_Spectrum.png`' )
            plt.savefig( 'figs/Wavelength_Calibrated_1-D_Spectrum.png', dpi = 144 )
        else:
            print( '' )
        if show:
            print( '[Wavelength calibration] Show plot' )
            plt.show()
        plt.close()

    if inverse: wav_caled = wav_caled[::-1]

    if not os.path.exists( 'cal' ): os.makedirs( 'cal' )
    file_path = F'cal/WaveSoln.dat'
    print( F'[Wavelength calibration] Write calibrated wavelength to `{file_path}`' )
    np.savetxt( file_path, np.vstack([ wav_caled, arc ]).T, fmt = '%15.8e' )

    return wav_caled

def sens_func( wav, cnt, exp, seeing, sw, airmass, extfile, stdfile, wave_range, order, show = 1, save = 0 ):
    '''
    '''

    # Speed of light
    c = 2.99792458e+10 # [cm/s]
    
    inverse = False
    if wav[1] < wav[0]: wav = wav[::-1]; cnt = cnt[::-1]; inverse = True
    
    # Slit loss correction
    # --------------------
    print( F'[Sensitivity function] Slit loss correction: FWHM = {seeing:.2f} [px], Slit Width = {sw:.2f} [px]' )
    sl = get_slit_loss( fwhm = seeing, sw = sw )
    cnt = cnt / ( 1 - sl )
    
    # Atmospheric extinction correction
    # ---------------------------------
    print( '[Sensitivity function] Extinction correction: Load extinction coefficients' )
    wav_ext, ext = np.loadtxt( extfile ).T
    print( '[Sensitivity function] Extinction correction: Dereddening' )
    ext = interp1d( wav_ext, ext, kind = 'quadratic' )( wav )
    ext_corr_factor = 10**( 0.4 * airmass * ext )
    cnt = cnt * ext_corr_factor
    
    # Convert to counts/A/s
    # ---------------------
    # 1. Bandpass
    dwav = np.abs( np.diff( wav ) )
    dwav = np.hstack([ dwav[0], dwav, dwav[-1] ])
    dwav = ( dwav[:-1] + dwav[1:] ) / 2
    # 2. Convert
    cnt = cnt / dwav / exp

    # Sensitivity function
    # --------------------
    if not os.path.exists( os.path.join( os.getcwd(), stdfile ) ):
        if not os.path.exists( os.path.join( os.path.split( os.path.realpath( __file__ ) )[0], 'lib/onedstds/', stdfile ) ):
            raise FileNotFoundError( 'No standard file found.' )
            sys.exit()
        else:
            stdfile = os.path.join( os.path.split( os.path.realpath( __file__ ) )[0], 'lib/onedstds/', stdfile )
    else:
        stdfile = os.path.join( os.getcwd(), stdfile )
    # 1. Load standard spectrum
    print( '[Sensitivity function] Load archived standard spectrum' )
    stdspec = np.loadtxt( stdfile, skiprows = 1 ).T
    if stdspec.shape[0] == 3:
        wav_mag, mag, bp = stdspec
    elif stdspec.shape[0] == 2:
        wav_mag, mag = stdspec
        dwav_mag = np.abs( np.diff( wav_mag ) )
        dwav_mag = np.hstack([ dwav_mag[0], dwav_mag, dwav_mag[-1] ])
        bp = ( dwav_mag[:-1] + dwav_mag[1:] ) / 2
    flx_mod = 10**( -0.4 * mag ) * 3631e-23 * c / wav_mag**2  * 1e8 # [erg/cm2/s/A]
    # 2. Comparision
    print( '[Sensitivity function] Comparision' )
    flx_obs = np.zeros( flx_mod.shape[0] )
    bins = np.vstack([ wav_mag - bp/2, wav_mag + bp/2 ]).T
    for i, bin in enumerate( bins ):
        idx = ( bin[0] < wav ) & ( wav < bin[1] )
        edges = interp1d( wav, cnt, 'linear', bounds_error = False, fill_value = np.nan )( bin )
        flx_obs[i] = np.trapz( np.hstack([ edges[0], cnt[idx], edges[1] ]), x = np.hstack([ bin[0], wav[idx], bin[1] ]) ) / ( bin[1] - bin[0] )
    sen = flx_obs / flx_mod
    # 3. Filter NaN
    idx = np.isnan( sen )
    wav_sen = wav_mag[~idx]; sen = 2.5 * np.log10( sen[~idx] )
    
    print( '[Sensitivity function] Fitting' )
    mask = ( wave_range[0] < wav_sen ) & ( wav_sen < wave_range[1] )
    for i in range( 5 ):
        knots = np.r_[ ( wav_sen[mask][0], ) * ( order + 1 ), ( wav_sen[mask][-1], ) * ( order + 1 ) ]
        spl = make_lsq_spline( wav_sen[mask], sen[mask], t = knots, k = order )
        mask = mask & ~sigma_clip( sen - spl( wav_sen ), sigma = 0.5, maxiters = 1, masked = True ).mask
    
    sen_fit = spl( wav )

    # Plot
    print( '[Sensitivity function] Plot fitted sensitivity function', end = '' )
    fig, ax = plt.subplots( 2, 1, figsize = (12, 12) )
    # Sensitivity function
    ax[0].plot( wav_sen[ mask], sen[ mask], '+', color = 'black',     ms = 10 )
    ax[0].plot( wav_sen[~mask], sen[~mask], '+', color = 'lightgrey', ms = 10 )
    ax[0].plot( wav, sen_fit, '-', c = 'red', lw = 2 )
    # Settings
    ax[0].grid( axis = 'both', color = '0.95', zorder = 0 )
    ax[0].set_xlim( wav.min(), wav.max() )
    ax[0].tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 18 )
    ax[0].set_ylabel( 'Sensitivity Function', fontsize = 22 )
    ax[0].set_title( 'Sensitivity Function', fontsize = 24 )

    ax[1].plot( wav_sen[ mask], ( sen - spl( wav_sen ) )[ mask], '+', color = 'black',     ms = 10, zorder = 1 )
    ax[1].plot( wav_sen[~mask], ( sen - spl( wav_sen ) )[~mask], '+', color = 'lightgrey', ms = 10, zorder = 0 )
    ax[1].axhline( y = 0, c = 'red', ls = '--' , zorder = 2 )
    # Settings
    ax[1].grid( axis = 'both', color = '0.95' )
    ax[1].set_xlim( wav.min(), wav.max() )
    ax[1].set_ylim( ( sen - spl( wav_sen ) )[mask].min() * 0.8, ( sen - spl( wav_sen ) )[mask].max() * 1.2 )
    ax[1].tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 18 )
    ax[1].set_xlabel( 'Wavelength [$\\mathrm{\\AA}$]', fontsize = 22 ); 
    ax[1].set_ylabel( 'Residuals', fontsize = 22 )
    fig.align_ylabels()
    
    if save:
        print( ' to `figs/Sensitivity_Function.png`' )
        plt.savefig( 'figs/Sensitivity_Function.png', dpi = 144 )
    else:
        print( '' )

    if show:
        print( '[Sensitivity function] Show plot' )
        plt.show()
    plt.close()
    
    if inverse: sen_fit = sen_fit[::-1]; wav = wav[::-1]; cnt = cnt[::-1]
    
    # Write to file
    if not os.path.exists( 'cal' ): os.makedirs( 'cal' )
    file_path = F'cal/SensFunc.dat'
    print( F'[Sensitivity function] Write sensitivity function to `{file_path}`' )
    np.savetxt( file_path, np.vstack([ wav, sen_fit ]).T, fmt = '%15.8e' )
    
    return sen_fit

def flux_calibration( wav, cnts, cnts_err, sens, exp, airmass, extfile, spectype, show = 1, save = 0 ):
    '''
    '''

    print( '[Flux calibration] Extinction correction: Load extinction coefficients' )
    wav_ext, ext = np.loadtxt( extfile ).T
    ext = interp1d( wav_ext, ext, kind = 'quadratic' )( wav )
    ext_corr_factor = 10**( 0.4 * airmass * ext )
    
    # Bandpass
    dwav = np.abs( np.diff( wav ) )
    dwav = np.hstack([ dwav[0], dwav, dwav[-1] ])
    dwav = ( dwav[:-1] + dwav[1:] ) / 2

    for i, ( cnt, cnt_err ) in enumerate( zip( cnts, cnts_err ) ):

        # Atmospheric extinction correction
        print( '[Flux calibration] Extinction correction: Dereddening' )
        cnt     = cnt     * ext_corr_factor
        cnt_err = cnt_err * ext_corr_factor

        # Convert to counts/A/s
        cnt     = cnt     / dwav / exp
        cnt_err = cnt_err / dwav / exp
        
        # Calibration
        flxs     = cnt     * 10**( -0.4 * sens )  # [erg/s/cm2/A]
        flxs_err = cnt_err * 10**( -0.4 * sens )  # [erg/s/cm2/A]

        # Plot
        print( F'[Flux calibration] Plot flux calibrated {i+1} of {len(cnts)} {spectype} spectrum', end = '' )
        fig, ax = plt.subplots( 1, 1, figsize = (12, 6) )
        # Spectrum
        for flx, flx_err in zip( flxs, flxs_err ):
            ax.step( wav, flx / 1e-14, 'b-', where = 'mid' )
            ax.fill_between( wav, ( flx - flx_err ) / 1e-14, ( flx + flx_err ) / 1e-14, color = 'b', alpha = 0.1 )
        # Settings
        ax.grid( axis = 'both', color = '0.95' )
        ax.set_xlim( wav.min(), wav.max() )
        ax.tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 18 )
        ax.set_xlabel( 'Wavelength [$\\mathrm{\\AA}$]', fontsize = 22 ); 
        ax.set_ylabel( 'Flux [$\\mathrm{10^{-14}\ erg\ s^{-1}\ cm^{-2}\ \AA^{-1}}$]', fontsize = 22 )
        ax.set_title( F'Flux Calibrated {spectype} Spectrum {i+1}', fontsize = 24 )
        
        if save:
            print( F' to `{ F"figs/{spectype}_FlxCal_{str(i+1).zfill(len(str(len(cnts))))}.png" }`' )
            plt.savefig( F'figs/{spectype}_FlxCal_{str(i+1).zfill(len(str(len(cnts))))}.png', dpi = 144 )
        else:
            print( '' )
        if show:
            plt.show()
        plt.close()

        # Write to file
        if not os.path.exists( 'cal' ): os.makedirs( 'cal' )
        file_path = F'cal/{spectype}_FlxCal_{str(i+1).zfill(len(str(len(cnts))))}.dat'
        print( F'[Flux calibration] Write flux calibrated {i+1} of {len(cnts)} {spectype} spectrum to `{file_path}`' )
        with open( file_path, 'w' ) as f:
            f.write( F'# Number of apertures: {flxs.shape[0]}\n' )
            f.write( F'# Flux unit: erg/s/cm2/A\n' )
            if wav[1] < wav[0]: 
                np.savetxt( f, np.vstack([ wav, flxs, flxs_err ]).T[::-1], fmt = '%15.8e' )
            else: 
                np.savetxt( f, np.vstack([ wav, flxs, flxs_err ]).T, fmt = '%15.8e' )