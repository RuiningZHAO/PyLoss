import os

try:
    from astropy.io import fits
    from astropy.time import Time
    from astropy.version import version
except ImportError:
    print( 'Module `astropy` not found. Please install with: pip install astropy' )
    sys.exit()

try:
    import numpy as np
except ImportError:
    print( 'Module `numpy` not found. Please install with: pip install numpy' )
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

def Login_info( version ):
    '''
    Login infomation
    '''
    from datetime import datetime
    
    year = datetime.now().strftime( '%Y' )

    print( '\n' )
    print( '< Pipeline for Long-slit Spectra (Pilloss) >'.center( 76 ) )
    print( F'Copyright 2018-{year} Ruining ZHAO'.center( 76 ) )
    print( F'Version {version}'.center( 76 ) )
    print( 'Release @ July 8, 2020'.center( 76 ) )
    print( '\n' )

def load_lists( listnames ):
    '''
    A function to load lists.

    Parameters
    ----------
    listnames : list
        Names of the lists.

    Returns
    -------
    lists : dict
        Loaded lists.
    '''
    
    lists = dict()
    for listname in listnames:
        with open( F'lists/{listname}.lst' ) as f:
            lists[listname] = [ line.strip() for line in f.readlines() ]
    
    return lists

def load_image( filename, extension ):
    '''
    '''
    hdr  = fits.getheader( filename, extension )
    data = fits.getdata( filename, extension, header = False )

    return data, hdr

def load_images( filelist, extension ):
    '''
    '''
    hdrlst = list()
    for i, file in enumerate( filelist ):
        img, hdr = load_image( os.path.join( os.getcwd(), file ), extension )
        hdrlst.append( hdr )
        if i == 0: imglst = np.zeros([ len( filelist ), img.shape[0], img.shape[1] ])
        imglst[i] = img
    errlst = np.sqrt( imglst )
    
    return imglst, errlst, hdrlst

def imcmb( imgs, errs, hdrs, method = 'median', reject = 'sigma_clip' ):
    '''
    '''
    
    # Rejection
    # ---------
    # 1. Min-max pre-rejection
    maxidx, minidx = imgs.argmax( axis = 0 ), imgs.argmin( axis = 0 )
    idx = np.tile( np.arange( imgs.shape[0] ), [imgs.shape[1], imgs.shape[2], 1] ).transpose( 2, 0, 1 )
    mask = ( idx == maxidx ) | ( idx == minidx )
    # 2. Rejection
    if reject == 'sigma_clip':
        # Median & standard deviation after min-max rejection
        med = np.ma.median( np.ma.array( imgs, mask = mask ), axis = 0 )
        sig = np.ma.std(    np.ma.array( imgs, mask = mask ), axis = 0, ddof = 1 )
        # Sigma clipping
        mask = np.abs( imgs - med ) > 3 * sig

    elif reject == 'minmax':
        pass

    else:
        raise ValueError( 'Unknown rejection. Only `minmax` and `sigma_clip` are available.' )

    # Combine
    # -------
    # Uncertainty
    err = ( np.ma.sqrt( np.ma.array( errs**2, mask = mask ).sum( axis = 0 ) ) / np.sum( ~mask, axis = 0 ) )
    err = err.filled( np.ma.median( err ) )
    # Combine
    if method == 'median':
        img = np.ma.median( np.ma.array( imgs, mask = mask ), axis = 0 )
    elif method == 'mean':
        img = np.ma.mean(   np.ma.array( imgs, mask = mask ), axis = 0 )
    else:
        raise ValueError( 'Unknown method. Only `median` and `mean` are available.' )
    img = img.filled( 99999 ) # ?????
    
    # Header
    hdrs[0]['BANDID1']  = F'Combined image - {reject.replace("_", " ")}, {method}'
    hdrs[0]['BANDID2']  = F'Uncertainty - {reject}'
    hdrs[0]['NCOMBINE'] = ( imgs.shape[0], 'Number of files combined' )
        
    return img, err, hdrs[0]
    
def Gaussian( x, A, mu, sig ):
    '''
    '''
    
    y = A * np.exp( -( x - mu )**2 / ( 2 * sig**2 ) )

    return y

def concatenate( specs, idx, mode ):
    '''
    '''
    
    print( F'[Concatenate] At the index of {idx} (mode = {mode})' )

    if mode == 1:
        pass
    elif mode == 2:
        specs = specs[::-1]

    spec = np.hstack([ specs[0][:idx], specs[1][idx:] ])

    return spec

# def polyfit( x, b, order ):
#     '''
#     '''
#     import numpy as np

#     X = np.array([ x**(order-i) for i in range( order+1 ) ] ).T
#     A = np.dot( np.linalg.inv( np.dot( X.T, X ) ), np.dot( X.T, b ) )

#     return A

# def polyval( x, A ):
#     '''
#     '''
#     import numpy as np

#     order = A.shape[0] - 1
#     X = np.array([ x**(order-i) for i in range( order+1 ) ] ).T

#     return np.dot( X, A )

def plot2d( img, title, slit_along, show = 1, save = 1 ):
    '''
    '''

    fig, ax = plt.subplots( 1, 1, figsize = ( 10, 8 ) )
    fig.subplots_adjust( right = 0.8 )

    # Image
    im = ax.imshow( img, cmap = 'Greys_r', origin = 'lower', extent = ( 0.5, img.shape[1]+0.5, 0.5, img.shape[0]+0.5) )
    # Colorbar
    cax = fig.add_axes([ ax.get_position().x1 + 0.02, ax.get_position().y0, 0.04, ax.get_position().height ])
    cb = fig.colorbar( im, cax = cax )
    # Settings
    ax.tick_params( which = 'major', direction = 'in', color = 'w', top = True, right = True, length = 7, width = 1.5, labelsize = 18 )
    if slit_along == 'col': ax.set_xlabel( 'Dispersion', fontsize = 22 ); ax.set_ylabel( 'Slit', fontsize = 22 )
    if slit_along == 'row': ax.set_xlabel( 'Slit', fontsize = 22 ); ax.set_ylabel( 'Dispersion', fontsize = 22 )
    ax.set_title( title, fontsize = 22 )
    cb.ax.tick_params( which = 'major', direction = 'in', color = 'w', right = True, length = 7, width = 1.5, labelsize = 18 )
    cb.ax.set_ylabel( 'ADU', fontsize = 22 )

    if show:
        print( '[Plotting] Show plot' )
        plt.show()

    if save:
        fig_path = 'figs'
        if not os.path.exists( fig_path ): os.makedirs( fig_path )
        print( F'[Plotting] Save to { os.path.join( fig_path, F"{title}.png".replace( " ", "_" ) ) }' )
        plt.savefig( os.path.join( fig_path, F'{title}.png'.replace( ' ', '_' ) ), dpi = 144 )

    plt.close()

    return None

def write( img, err, hdr, fit_path, title ):
    '''
    '''

    if not os.path.exists( fit_path ): os.makedirs( fit_path )

    print( F'[Write] Save to { os.path.join( fit_path, F"{title}.fits".replace( " ", "_" ) ) }' )
    # Edit header
    hdr['EXTEND'] = ( False, 'File may contain extensions' )
    hdr['ORIGIN'] = ( F'Astropy v{version}', 'FITS file originator' )
    hdr['DATE']   = ( Time.now().fits, 'Date FITS file was generated' )
    for key in ['BSCALE', 'BZERO']:
        if key in hdr.keys(): del hdr[key]
    # Write
    fits.writeto( os.path.join( fit_path, F"{title}.fits".replace( " ", "_" ) ), 
                  data      = [img.astype( np.float32 ), err.astype( np.float32 )], 
                  header    = hdr, 
                  overwrite = True )

# def plot3d(data, path, name, show = 1, save = 0):

#     import os
#     import numpy as np
#     import matplotlib.pyplot as plt
#     from mpl_toolkits.mplot3d import axes3d

#     figure_path = os.path.join( path, 'figures' )

#     if not os.path.exists( figure_path ):
#         os.system( F'mkdir {figure_path}' )

#     fig = plt.figure( figsize = (10, 6) )
#     ax = fig.add_subplot( 111, projection = '3d' )
#     X, Y = np.meshgrid( np.arange( len(data[0]) ), np.arange( len(data) ) )
#     ax.plot_wireframe( X, Y, data, color = 'k', lw = 0.8, rstride = 50, cstride = 50, alpha = 0.5 )
#     ax.set_xlabel( 'Column', labelpad = 20, fontsize = 16 )
#     ax.set_ylabel( 'Row',    labelpad = 20, fontsize = 16 )
#     ax.set_zlabel( 'ADU',    labelpad = 20, fontsize = 16 )

#     if save:
#         plt.savefig( os.path.join( figure_path, F'{name}.png' ), dpi = 300 )

#     if show:
#         plt.show()

#     plt.close()

#     return None