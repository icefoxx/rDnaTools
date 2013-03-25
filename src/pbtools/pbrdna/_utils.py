import os
import sys

def fileExists( filename ):
    return os.path.exists(filename) and (os.path.getsize(filename) > 0)

def isExe( filePath ):
    if filePath is None:
        return False
    return os.path.isfile(filePath) and os.access(filePath, os.X_OK)

def getZmw( read ):
    parts = read.split('/')
    return '/'.join(parts[0:2])

def returnEmpty():
    return []

def createDirectory( dirName ):
    try:
        os.mkdir( dirName )
    except OSError:
        pass
    if not os.path.isdir( dirName ):
        raise OSError('Could not create directory "%s"!' % dirName)
    return dirName

def validateInputFile( fileName, allowedSuffixes ):
    # First we check whether the input file has a valid suffix
    try: 
        assert any( [fileName.endswith(suffix) for suffix in allowedSuffixes] )
    except AssertionError:
        raise ValueError('File does not have an allowed suffix! %s' % \
                                                        allowedSuffixes)
    # Next we check whether the input file exists where specified
    try: 
        assert fileExists( fileName )
    except AssertionError:
        raise OSError('Input file does not exist!')
    # Finally we return the absolute path to the file
    return os.path.abspath( fileName )

def validateOutputFile( fileName ):
    if fileName in [sys.stdout, sys.stderr]:
        return fileName
    return os.path.abspath( fileName )

def which(program):
    """
    Find and return path to local executables  
    """
    fpath, fname = os.path.split(program)
    if fpath:
        if isExe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exeFile = os.path.join(path, program)
            if isExe(exeFile):
                return exeFile
    return None
