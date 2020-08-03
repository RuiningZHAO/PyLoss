import os, sys
from copy import deepcopy

try:
    import numpy as np
except ImportError:
    print( 'Module `numpy` not found. Please install with: pip install numpy' )
    sys.exit()

try:
    from astropy.time import Time
    from astropy.stats import sigma_clip
except ImportError:
    print( 'Module `astropy` not found. Please install with: pip install astropy' )
    sys.exit()

try:
    from scipy.signal import find_peaks
    from scipy.optimize import curve_fit
    from scipy.ndimage import gaussian_filter
    from scipy.interpolate import interp1d, griddata, make_lsq_spline
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
    from lacosmic import lacosmic
except ImportError:
    print( 'Module `lacosmic` not found. Please install with: pip install lacosmic' )
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

from .utils import load_images, Gaussian, imcmb, plot2d, write

class CCDDataList( object ):
    '''
    '''
    def __init__( self, filelist, extension, slit_along, rdnoise, gain ):
        '''
        A method

        Parameters
        ----------
        rdnoise : float
            Read out noise.
        gain : float
            Photon gain.
        '''
        self.filelist   = filelist
        self.extension  = extension
        self.slit_along = slit_along
        self.rdnoise    = rdnoise
        self.gain       = gain
        self.num        = len( filelist )
        self.combined   = False
        
        # Load images
        print( F'[Loading images] {self.num} images in total' )
        self.imgs, self.errs, self.hdrs = load_images( self.filelist, self.extension )
        
        return None

    def trim( self, col_range, row_range ):
        '''
        A method to trim images.

        Parameters
        ----------
        col_range : array_like
            Column range.
        row_range : array_like
            Row range.
        '''
        
        # Trim
        print( F'[Trimming] Trim section: column range: {col_range}, row range: {row_range}' )
        self.imgs = self.imgs[:, (row_range[0]-1):row_range[1], (col_range[0]-1):col_range[1]]
        self.errs = self.errs[:, (row_range[0]-1):row_range[1], (col_range[0]-1):col_range[1]]
        
        # Edit header
        for hdr in self.hdrs:
            hdr['WCSDIM']   = 2
            hdr['LTM1_1']   = 1.
            hdr['LTM2_2']   = 1.
            hdr['WAT0_001'] = 'system=physical'
            hdr['WAT1_001'] = 'wtype=linear'
            hdr['WAT2_001'] = 'wtype=linear'
            hdr['TRIM']     = F'{Time.now().fits} Trim data section is [{col_range[0]}:{col_range[1]},{row_range[0]}:{row_range[1]}]'
            hdr['CCDSEC']   = F'[{col_range[0]}:{col_range[1]},{row_range[0]}:{row_range[1]}]'
            hdr['BIASSEC']  = F'[1:{col_range[1]-col_range[0]+1},1:{row_range[1]-row_range[0]+1}]'
            hdr['LTV1']     = 1 - col_range[0]
            hdr['LTV2']     = 1 - row_range[0]

        return deepcopy( self.imgs ), deepcopy( self.errs ), deepcopy( self.hdrs )

    def combine( self, method, reject ):
        '''
        A method to combine images.
        
        Parameters
        ----------
        method : string
        reject : string
        '''

        # Combine
        print( F'[Combining] #images: {self.num}, method: {method}, reject: {reject}' )
        self.img, self.err, self.hdr = imcmb( self.imgs, self.errs, self.hdrs, method, reject )
        self.combined = True

        # Edit header
        self.hdr['COMMENT'] = 'Combined from '+ ', '.join( self.filelist )
        
        return deepcopy( self.img ), deepcopy( self.err ), deepcopy( self.hdr )
        
    def info( self, combined = False ):
        '''
        A method to print statistics.

        Parameters
        ----------
        combined : bool
            If True, statistics for combined image will be printed.
        '''

        tab = PrettyTable( hrules = HEADER, vrules = NONE )
        if not combined:
            tab.field_names = [ 'FILE NAME', 'SHAPE', 'MEAN', 'STDDEV', 'MIN', 'MAX' ]
            for fname, img, err in zip( self.filelist, self.imgs, self.errs ):
                tab.add_row([ fname,         img.shape, round( img.mean(), 2 ), round( img.std(ddof=1), 2 ), int( img.min() ), int( img.max() ) ])
                tab.add_row([ 'Uncertainty',      '\"', round( err.mean(), 2 ), round( err.std(ddof=1), 2 ), round( err.min(), 2 ), round( err.max(), 2 ) ])

        elif self.combined:
            tab.field_names = [ 'SHAPE', 'MEAN', 'STDDEV', 'MIN', 'MAX' ]
            tab.add_row([ self.img.shape, round( self.img.mean(), 2 ), round( self.img.std(ddof=1), 2 ), round( self.img.min(), 2 ), round( self.img.max(), 2 ) ])
            tab.add_row([           '\"', round( self.err.mean(), 2 ), round( self.err.std(ddof=1), 2 ), round( self.err.min(), 2 ), round( self.err.max(), 2 ) ])
        
        else:
            raise RuntimeError( 'The method `combine` should be called first if `combined` is True.' )
            sys.exit()

        for col in ['MEAN', 'STDDEV', 'MIN', 'MAX']: tab.align[col] = 'r'
        print( '\n' + tab.get_string() + '\n' )
    
    def plot( self, title, combined = False, show = True, save = False ):
        '''
        '''

        if not combined:
            for i, img in enumerate( self.imgs ):
                plot2d( img, F'{ title } { str(i+1).zfill(len(str(self.num))) }', self.slit_along, show, save )
        elif self.combined:
            plot2d( self.img, title, self.slit_along, vmin = None, show = show, save = save )
        else:
            raise RuntimeError( 'The method `combine` should be called first if `combined` is True.' )
            sys.exit()

    def write( self, path, title, combined = False ):
        '''
        '''

        if not combined:
            for i, ( img, err, hdr ) in enumerate( zip( self.imgs, self.errs, self.hdrs ) ):
                write( img, err, hdr, path, F'{ title } { str(i+1).zfill(len(str(self.num)))}' )

        elif self.combined:
            write( self.img, self.err, self.hdr, path, title )

        else:
            raise RuntimeError( 'The method `combine` should be called first if `combined` is True.' )
            sys.exit()

class Corrector( object ):
    '''
    '''
    def __init__( self, imgs, errs, hdrs, slit_along, rdnoise, gain ):
        '''
        '''
        self.imgs       = deepcopy( imgs )
        self.errs       = deepcopy( errs )
        self.hdrs       = deepcopy( hdrs )
        self.slit_along = slit_along
        self.rdnoise    = rdnoise
        self.gain       = gain
        self.num        = self.imgs.shape[0]
        
        return None

    def bias_subtraction( self, bias, berr ):
        '''
        '''

        print( '[Image correction] Bias subtraction' )
        # Bias correction
        self.imgs -= bias
        # Uncertainty
        self.errs = np.sqrt( self.errs**2 + berr**2 )

        # Header
        for hdr in self.hdrs:
            hdr['ZEROCOR'] = ( round( bias.mean(), 0 ), 'Zero level corrected' )

        return deepcopy( self.imgs ), deepcopy( self.errs ), deepcopy( self.hdrs )

    def normalization( self, sigma ):
        '''
        A method to normalize flatfield.

        Parameters
        ----------
        sigma : array_like
            sigma used in smoothing.
        '''

        if self.num > 1: 
            raise RuntimeError( 'This method is designed for flat-field normalization. Only one master flat-field array should be passed in.' )
            sys.exit()

        print( '[Image correction] Flat normalization' )
        # Smoothed Flat
        # -------------
        # 1. Initial smooth
        img_smth = gaussian_filter( self.imgs[0], sigma = sigma, order = 0, mode = 'mirror' )
        err_smth = gaussian_filter( self.errs[0], sigma = sigma, order = 0, mode = 'mirror' )
        # 2. 2-D Interpolate
        img_norm = self.imgs[0] / img_smth
        mask = np.abs( img_norm - 1 ) > 3 * np.std( img_norm - 1, ddof = 1 )
        x = np.arange( self.imgs[0].shape[1] ); y = np.arange( self.imgs[0].shape[0] )
        X, Y = np.meshgrid( x, y )
        img_intp = griddata( ( Y[~mask], X[~mask] ), self.imgs[0][~mask], ( Y, X ), method = 'nearest' )
        # 3. Second smooth
        img_smth = gaussian_filter( img_intp, sigma = sigma, order = 0, mode = 'mirror' )

        plot2d( mask.astype(int), 'Mask', self.slit_along, show = 1, save = 1 )

        # Edit header
        self.hdrs[0]['BANDID1'] = F'Normalized flat - 2-D gaussian_filter, sigma = { sigma }'
        self.hdrs[0]['BANDID2'] = 'Uncertainty'

        # Uncertainty
        self.errs[0] = np.sqrt( ( self.errs[0] / img_smth )**2 + ( self.imgs[0] * err_smth / img_smth**2 )**2 )
        # Normalization
        self.imgs[0] = self.imgs[0] / img_smth

        return ( deepcopy( self.imgs[0] ), deepcopy( self.errs[0] ) ), ( deepcopy( img_smth ), deepcopy( err_smth ) )

    def flat_correction( self, flat, ferr ):
        '''
        '''

        print( '[Image correction] Flat-fielding' )
        # Uncertainty
        self.errs = np.sqrt( ( self.errs / flat )**2 + ( self.imgs * ferr / flat**2 )**2 )
        # Flat correction
        self.imgs = self.imgs / flat

        # Header
        for hdr in self.hdrs:
            hdr['BANDID1'] = 'Preprocessed image - bias subtracted, flat-field corrected'
            hdr['BANDID2'] = 'Uncertainty'
            hdr['FLATCOR'] = ( round( flat.mean(), 0 ), 'Flat-field corrected' )
        
        return deepcopy( self.imgs ), deepcopy( self.errs ), deepcopy( self.hdrs )

    def maskcr( self, method ):
        '''
        A method to mask cosmic ray pixels.

        Parameters
        ----------
        method : `Laplacian` or `median`
            If len(imglist) > 2, two methods are available, which are 
            Laplacian Edge Detection and median with sigma clipping
            denoted as `Laplacian` and `median` respectively. If 
            len(imglist) <= 2, Laplacian Edge Detection is used.
        '''

        if self.num <= 2: method == 'Laplacian'

        if method == 'Laplacian':
            for i in range( self.num ):
                print( F'[Removing bad pixels] Process {i+1} of { self.num } images using { method } method [Cosmic ray & hot pixels]\n' )
                self.imgs[i], _ = lacosmic( self.imgs[i], 
                                            contrast           = 1, 
                                            cr_threshold       = 4.5, 
                                            neighbor_threshold = 0.5, 
                                            readnoise          = self.rdnoise, 
                                            effective_gain     = self.gain )
#                 if dead_px:
#                     print( F'\n[Removing bad pixels] Process {i+1} of { self.num } images using { method } method [Dead pixels]\n' )
#                     img = self.imgs[i]
#                     if self.slit_along == 'row': img = img.T
#                     # Mask bad pixels
#                     img_smth = gaussian_filter( img, ( 100, 0 ) )
#                     img_resi = img - img_smth
#                     mask = img_resi < -4 * np.std( img_resi, axis = 0, ddof = 1 )[np.newaxis, :]
#                     # 2-D Interpolate
#                     img = np.ma.array( img, mask = mask )
#                     x = np.arange( img.shape[1] ); y = np.arange( img.shape[0] )
#                     X, Y = np.meshgrid( x, y )
#                     img_itp = griddata( ( Y[~img.mask], X[~img.mask] ), img[~img.mask], ( Y, X ), method = 'linear' )
#                     if self.slit_along == 'row': img_itp = img_itp.T
#                     self.imgs[i] = img_itp

                self.hdrs[i]['COMMENT'] = 'Bad pixels removed using Laplacian Edge Detection'
                print( '' )

#         elif method == 'median':
#             print( F'[Remove Cosmic Rays] Process { self.num } { imagetype } images using { method } method\n' )
#             data = np.array([ img for img in self.imgs ])
#             mask = np.zeros( data.shape, dtype = bool )
#             for i in range( 3 ):
#                 # Median & standard deviation after maximum rejection
#                 maximum = data.max( axis = 0 )
#                 mask = mask | ( data == maximum )
#                 median = np.ma.median( np.ma.array( data, mask = mask ), axis = 0 )
#                 std    = np.ma.std(    np.ma.array( data, mask = mask ), axis = 0, ddof = 1 )
#                 # Sigma clipping & replace cosmic ray pixels with median value
#                 mask = ( data - median ) > 3 * std
#                 data_median = np.ma.median( np.ma.array( data, mask = mask ), axis = 0 )
#                 data[mask] = np.tile( data_median, [ data.shape[0], 1, 1 ] )[mask]
#             for i in range( len(imglist) ):
#                 print( F'{i+1} of {len(imglist)} images: {mask[i].sum()} pixels masked [median with sigma clipping]\n' )
#                 imglist[i][0] = data[i]
#                 hdrlist[i]['COMMENT'] = 'Cosmic ray pixels removed using median with sigma clipping'

        else:
            raise ValueError( F'Unknown method `{method}`. Only `Laplacian` are available.' )

    def get_zp_shift( self, order, show, save ):
        '''
        A method to get zeropoint shift along dispersion axis.

        Parameters
        ----------
        order : int
            Order used in polyfit.
        show : bool
            If `True`, the plot will be shown.
        save : bool
            If `True`, the plot will be written to file.

        Returns
        -------
        self.shift : array_like
            Fitted zeropoint shift.
        '''

        # Derive zeropoint shift
        # ----------------------
        shift = np.array([])
        for k, img in enumerate( self.imgs ):

            if self.slit_along == 'col': data = img * 1
            if self.slit_along == 'row': data = img.T

            # Find peaks
            # ----------
            print( F'[Zeropoint correction] Seeking for peaks in the {k+1}/{self.num} images' )
            peaks, properties = find_peaks( data.mean( axis = 0 ) / data.mean( axis = 0 ).max(), 
                                            height = ( 0.3, 0.8 ), 
                                            distance = data.shape[1] // 10,
                                            width = 0 )
            # Print peak info.
            tab = PrettyTable( hrules = HEADER, vrules = NONE )
            tab.field_names = [ 'CENTER', 'HEIGHT', 'WIDTH' ]
            for i, peak in enumerate( peaks ):
                tab.add_row([ peak, round( properties['peak_heights'][i], 2 ), int( properties['widths'][i] ) ])
            print( '\n' + tab.get_string() + '\n' )

            # Gaussian fitting
            x = np.arange( data.shape[1] )
            mu = np.zeros([ data.shape[0], len( peaks ) ])
            for i in range( data.shape[0] ):
                print( F'\r[Zeropoint correction] Peak center fitting in the ({i+1}/{data.shape[0]})-th row of {k+1}/{self.imgs.shape[0]} images', end = '', flush = True )
                for j, peak in enumerate( peaks ):
                    idxmin, idxmax = peak - int( properties['widths'][j] ), peak + int( properties['widths'][j] )
                    popt, pcov = curve_fit( Gaussian, 
                                            x[idxmin:idxmax], 
                                            data[i, idxmin:idxmax], 
                                            bounds = ( [data[i, idxmin:idxmax].max()*0.5, x[idxmin:idxmax].min(), 0 ], 
                                                       [data[i, idxmin:idxmax].max()*1.5, x[idxmin:idxmax].max(), x[idxmin:idxmax].shape[0] ] ) )
                    mu[i, j] = popt[1]
            # Convert to relative shift
            shift = np.hstack([ shift.reshape( mu.shape[0], -1 ), ( mu - mu[mu.shape[0]//2] ) ])
            print( '' )
        # Mean
        shift_mean = shift.mean( axis = 1 )

        # Polyfit to shift curve
        # ----------------------
        print( F'[Zeropoint correction] Fitting zeropoint shift iteratively (order = {order}, 5 iterations, nsigma = 1)' )
        y = np.arange( shift_mean.shape[0] ) + 1
        mask = np.ones( shift_mean.shape[0], dtype = bool )
        for k in range( 5 ):
            knots = np.r_[ ( y[mask][0], ) * ( order + 1 ), ( y[mask][-1], ) * ( order + 1 ) ]
            spl = make_lsq_spline( y[mask], shift_mean[mask], t = knots, k = order )
            mask = mask & ~sigma_clip( shift_mean - spl( y ), sigma = 1, maxiters = 1, masked = True ).mask
#             p = np.poly1d( np.polyfit( y[mask], shift_mean[mask], order ) )
#             shift_fit = p( y )
#             mask = mask & ( ~sigma_clip( shift_mean - shift_fit, sigma = 1, maxiters = 1, masked = True ).mask )
        self.shift = spl( y )

        # Write to file
        # -------------
        if not os.path.exists( 'bak' ): os.makedirs( 'bak' )
        np.savetxt( 'bak/zp_shift.dat', self.shift, fmt = '%15.8e' )

        # Plot
        # ----
        fig, ax = plt.subplots( 2, 1, figsize = ( 10, 8 ) )
        fig.subplots_adjust( hspace = 0 )
        
        # Shifts
        for i in range( shift.shape[1] ):
            ax[0].plot( y, shift[:, i], 'r+' )
        ax[0].plot( y[~mask], shift_mean[~mask], '+', c = 'grey' )
        ax[0].plot( y[mask], shift_mean[mask], 'k+' )
        # Fitted shifts
        ax[0].plot( self.shift, '-', c = 'yellow', lw = 2 )
        # Settings
        ax[0].set_xlim( y.min(), y.max() )
        ax[0].tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 18 )
        ax[0].set_xticklabels([])
        ax[0].set_ylabel( 'Displacement [px]', fontsize = 22 )
        ax[0].set_title( 'Zeropoint Shift Curve Fitting', fontsize = 24 )
        
        # Residuals
        ax[1].plot( y[mask],  shift_mean[mask]  - self.shift[mask], 'k+' )
        ax[1].plot( y[~mask], shift_mean[~mask] - self.shift[~mask], '+', c = 'grey' )
        # Settings
        ax[1].axhline( y = 0, ls = '--', c = 'yellow', lw = 2 )
        ax[1].set_xlim( y.min(), y.max() )
        ax[1].tick_params( which = 'major', direction = 'in', top = True, right = True, length = 7, width = 1.5, labelsize = 18 )
        ax[1].set_xlabel( 'Slit', fontsize = 22 )
        ax[1].set_ylabel( 'Residuals [px]', fontsize = 22 )
        fig.align_ylabels()
        
        if save:
            fig_path = 'figs'
            print( F'[Zeropoint correction] Plot zeropoint shift curve fitting to `{ os.path.join( fig_path, "Zeropoint_shift_fitting.png" ) }`' )
            plt.savefig( os.path.join( fig_path, 'Zeropoint_shift_curve_fitting.png' ), dpi = 144 )
        if show:
            print( F'[Zeropoint correction] Show plot' )
            plt.show()
        plt.close()

        return deepcopy( self.shift )

    def zp_shift_correction( self, shift ):
        '''
        '''

        for i, ( img, err, hdr ) in enumerate( zip( self.imgs, self.errs, self.hdrs ) ):

            print( F'[Zeropoint correction] Correction for {i+1} of { self.num } images' )
            if self.slit_along == 'row': data = img.T; erro = err.T
            if self.slit_along == 'col': data = img*1; erro = err*1

            # Resample
            # --------
            new_data = np.zeros( data.shape )
            new_erro = np.zeros( erro.shape )
            x = np.arange( data.shape[1] )
            for j in range( data.shape[0] ):
                print( F'\r[Zeropoint correction] Resampling ({j+1}/{data.shape[0]})-th row', end = '', flush = True )
                if shift[j] <= 0: fill = data[j][ 0]
                if shift[j] >  0: fill = data[j][-1]
                new_data[j], new_erro[j] = spectres( x, x - shift[j], data[j], erro[j], fill, verbose = False )
            print( '' )

            if self.slit_along == 'row': new_data = new_data.T; new_erro = new_erro.T
            
            self.imgs[i], self.errs[i] = new_data, new_erro
            self.hdrs[i]['COMMENT'] = 'Zeropoint shifted'
        
        return deepcopy( self.imgs ), deepcopy( self.errs ), deepcopy( self.hdrs )

    def plot( self, title, show = True, save = False ):
        '''
        '''

        for i, img in enumerate( self.imgs ):
            if self.num > 1:
                ttl = F'{ title } { str(i+1).zfill(len(str(self.num))) }'
            elif self.num == 1:
                ttl = title
            vmin = np.median( img ) - 3 * np.std( img - np.median( img ), ddof = 1 )
            plot2d( img, ttl, self.slit_along, vmin, show, save )

    def write( self, path, title ):
        '''
        '''

        for i, ( img, err, hdr ) in enumerate( zip( self.imgs, self.errs, self.hdrs ) ):
            if self.num > 1:
                ttl = F'{ title } { str(i+1).zfill(len(str(self.num))) }'
            elif self.num == 1:
                ttl = title
            write( img, err, hdr, path, ttl )

# if __name__ == '__main__':
#     pass
