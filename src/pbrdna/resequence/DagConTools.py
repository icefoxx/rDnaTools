#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this 
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors 
#   may be used to endorse or promote products derived from this software 
#   without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$

import os
import sys
import logging
import subprocess

from pbrdna._utils import createDirectory, validateInputFile, validateOutputFile
from pbrdna._utils import which

SCRIPT_CHOICES = ['gcon.py']
MODE_CHOICES = ['r', 'd']

class DagConRunner(object):
    """
    A tool for resequencing clusters of rDNA sequences with 
    the Gcon.py script from PB-DagCon
    """

    ##########################
    # Initialization Methods #
    ##########################

    def __init__(self, script, mode=None):
        self.script = script
        self.mode   = mode
        # Validate the settings
        self.validateSettings()

    def validateSettings(self):
        # Check the values of the specified script
        try:
            assert self.script in SCRIPT_CHOICES
        except:
            raise ValueError("Script is not a recognized DagCon script!")
        # If it's gcon, check that a valid mode was specified
        if self.script == 'gcon.py':
            try:
                assert self.mode in MODE_CHOICES
            except:
                raise ValueError("Gcon.py runners must specify mode 'r' or 'd'")
        # Finally, if the script and options pass, find the absolute pathes
        self.executable = which( self.script )

    ####################
    # Instance Methods #
    ####################

    def runGcon(self, inputFile, outputFile, refFile=None):
        if outputFile is None:
            outputFile = self.getOutputFileName( inputFile )
        if self.mode == 'r':
            assert refFile is not None
            p = subprocess.Popen( [self.executable, 
                                   self.mode,
                                   inputFile,
                                   refFile,
                                   '-o', outputFile] )
            p.wait()
        elif self.mode == 'd':
            p = subprocess.Popen( [self.executable, 
                                   self.mode,
                                   inputFile,
                                   '-o', outputFile] )
            p.wait()
        return outputFile

    def getOutputFile(self, inputFile):
        path, filename = os.path.split( inputFile )
        root, ext = os.path.splitext( filename )
        outputFile = root + '_consensus.fa'
        return os.path.join( path, outputFile )

    def __call__(self, inputFile, refFile=None):
        outputFile = self.getOutputFile( inputFile )
        if os.path.exists( outputFile ):
            return outputFile
        elif self.script == 'gcon.py':
            return self.runGcon( inputFile, outputFile, refFile )
