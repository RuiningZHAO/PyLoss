import os
import errno
from shutil import copy2

def copy_input():
    '''
    '''

    copy2( os.path.join( os.path.split( os.path.realpath( __file__ ) )[0], 'lib/input.py' ), os.getcwd() )

def check_file( filename ):
    '''
    '''

    if os.path.exists( os.path.join( os.getcwd(), filename ) ):
        pass
    else:
        raise FileNotFoundError( errno.ENOENT, 
                                 os.strerror( errno.ENOENT ), 
                                 os.path.join( os.getcwd(), filename ) )

def check_folders( folders ):
    '''
    '''

    for folder in folders:
        path = os.path.join( os.getcwd(), folder )
        if not os.path.exists( path ): os.makedirs( path )
    
def main( listnames, folders ):

    copy_input()

    for listname in listnames:
        check_file( F'lists/{listname}.lst' )

    check_folders( folders )
    
if __name__ == '__main__':

    main( listnames = ['bias', 'flat', 'arc', 'obj', 'std'],
          folders = ['figs', 'corr', 'bak', 'output'] )