#!/usr/bin/env python

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

import sys
import logging

from pbcore.io import BasH5Reader
from FastqIO import FastqRecord, FastqWriter 

class BasH5Extractor(object):
    """
    A tool for extracting sequence data from a BasH5 file    
    """

    ##########################
    # Initialization Methods #
    ##########################

    def __init__(self, fileName=None, outputFile=None):
        if fileName is None:
            self.initializeFromArgs()
        else:
            self.initializeFromCall(fileName, outputFile)
        self.validateSettings()
        self.initializeReaders()

    def initializeFromArgs(self):
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('inputFile', metavar='FILE',
                            help="BasH5 or FOFN to extract from")
        parser.add_argument('-o', '--output', default=sys.stdout,
                            help="Specify a file to output the data to")
        group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument('--subreads',  action='store_true',
                            help="Output sequences from individual subreads")
        group.add_argument('--CCS', action='store_true',
                            help="Output sequences from CCS reads")
        args = parser.parse_args()
        self.__dict__.update( vars(args) )

    def initializeFromCall(self, fileName, outputFile):
        self.inputFile = fileName
        if outputFile is None:
            self.output = sys.stdout
        else:
            self.output = outputFile

    def validateSettings(self):
        try:
            assert self.inputFile.endswith('.bas.h5') or \
                   self.inputFile.endswith('.fofn')
        except:
            raise ValueError('Input files must be either FOFN or BasH5!')
        logging.info('Creating a BasH5Extractor for "%s"' % self.inputFile)
        logging.info('Outputing extracted reads to "%s"' % self.output)

    def initializeReaders(self):
        if self.inputFile.endswith('.bas.h5'):
            logging.info('Creating a BasH5Reader for "%s"' % self.inputFile)
            self.basH5Readers = [BasH5Reader( self.inputFile )]
        elif self.inputFile.endswith('.fofn'):
            self.basH5Readers = self.parseFofnFile( self.inputFile )

    def parseFofnFile(self, fofnFile):
        basH5Readers = []
        with open(fofnFile, 'r') as handle:
            for line in handle:
                fofnEntry = line.strip()
                try:
                    assert fofnEntry.endswith('.bas.h5')
                except:
                    raise ValueError('FOFN must contain only BasH5 files!')
                logging.info('Creating a BasH5Reader for "%s"' % fofnEntry)
                reader = BasH5Reader( fofnEntry )
                basH5Readers.append( reader )
        return basH5Readers

    #################
    # Class Methods #
    #################

    @classmethod
    def writeCcsFastq(cls, basH5Reader, fastqWriter):
        logging.info('Writing Fastq CCS reads from "%s"...' % basH5Reader.movieName)
        for zmw in basH5Reader:
            ccsRead = zmw.ccsRead()
            if ccsRead is not None:
                fastqRecord = FastqRecord(ccsRead.readName,
                                          ccsRead.basecalls(),
                                          ccsRead.qv('QualityValue'))
                fastqWriter.writeRecord( fastqRecord )

    @classmethod
    def writeSubreadFastq(cls, basH5Reader, fastqWriter):
        logging.info('Writing Fastq subreads from "%s"...' % basH5Reader.movieName)
        for zmw in basH5Reader:
            for subread in zmw.subreads():
                fastqRecord = FastqRecord(subread.readName,
                                          subread.basecalls(),
                                          subread.qv('QualityValue'))
                fastqWriter.writeRecord( fastqRecord )

    ####################
    # Instance Methods #
    ####################

    def outputCcsFastq(self):
        logging.info('Parsing Fastq CCS reads from input BAS.H5 files')
        with FastqWriter(self.output) as writer:
            for reader in self.basH5Readers:
                self.writeCcsFastq( reader, writer )

    def outputSubreadFastq(self):
        logging.info('Parsing Fastq subreads from input BAS.H5 files')
        with FastqWriter(self.output) as writer:
            for reader in self.basH5Readers:
                self.writeSubreadFastq( reader, writer )

    def __call__(self):
        if self.CCS:
            self.outputCcsiFastq()
        elif self.subreads:
            self.outputSubreadFastq()

if __name__ == '__main__':
    extractor = BasH5Extractor()
    extractor()
