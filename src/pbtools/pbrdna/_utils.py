import os
import sys

from collections import namedtuple

BlasrM1 = namedtuple('BlasrM1', ['qname', 'tname', 'qstrand', 'tstrand',
                                 'score', 'pctsimilarity', 
                                 'tstart', 'tend', 'tlength',
                                 'qstart', 'qend', 'qlength',
                                 'ncells'])

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

def validateInputFile( fileName, allowedSuffixes ):
    allowedSuffixes = [allowedSuffixes] if isinstance(str, type(allowedSuffixes)) \
                                        else allowedSuffixes
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

def validateExecutable( executable ):
    exePath = which( executable )
    try: 
        assert exePath is not None
    except AssertionError:
        raise ValueError('"%s" is not a valid executable!' % executable)
    return exePath

def validateInt( integer, minValue=None, maxValue=None ):
    try: # First we check whether the supplied parameter is an Int
        assert isinstance(integer, int)
    except AssertionError:
        raise TypeError('Parameter is not an Integer!')
    # If a minimum value is supplied, compare it to the Integer
    if minValue is not None:
        try:
            assert integer >= minValue
        except AssertionError:
            raise ValueError("Integer is less than Minimum Value!")
    # If a maximum value is supplied, compare it to the Integer
    if maxValue is not None:
        try:
            assert integer <= maxValue
        except AssertionError:
            raise ValueError("Integer is greater than Maximum Value!")
    return integer

def validateFloat( floating_point, minValue=None, maxValue=None ):
    try: # First we check whether the supplied parameter is an Int
        assert isinstance(floating_point, float)
    except AssertionError:
        raise TypeError('Parameter is not a Floating Point!')
    # If a minimum value is supplied, compare it to the Integer
    if minValue is not None:
        try:
            assert floating_point >= minValue
        except AssertionError:
            raise ValueError("Float is less than Minimum Value!")
    # If a maximum value is supplied, compare it to the Integer
    if maxValue is not None:
        try:
            assert floating_point <= maxValue
        except AssertionError:
            raise ValueError("Float is greater than Maximum Value!")
    return floating_point
