#!/usr/bin/env python

import sys
import os
import re
import logging

from numpy import where

from pbrdna.io.FastqIO import FastqReader, FastqRecord, FastqWriter

__version__ = "0.1"

MIN_QV = None
MIN_LENGTH = 100

class QualityTrimmer(object):
    """
    Tool for trimming low-quality bases from the ends of FASTQ files 
    """
    ##########################
    # Initialization Methods #
    ##########################

    def __init__(self, fastqFile=None, 
                       outputFile=None, 
                       minQv=None, 
                       minLength=None):
        if fastqFile is None:
            self.initializeFromArgs()
        else:
            self.initializeFromCall(fastqFile, outputFile, minQv, minLength)
        self.validateSettings()

    def initializeFromArgs(self):
        import argparse
        desc = 'A tool for trimming low-quality ends from FASTQ quality values'
        parser = argparse.ArgumentParser(description=desc)
        parser.add_argument('fastq', metavar='FASTQ_FILE',
                            help="File of FASTQ sequence data to align")
        parser.add_argument('-q', '--minimum_qv', metavar='INT',
                            type=int, default=MIN_QV, dest='minQv',
                            help="Minimum QV to leave un-masked")
        parser.add_argument('-l', '--minimum_length', metavar='INT',
                            type=int, default=MIN_LENGTH, dest='minLength',
                            help="Minimum length to require " + \
                                 "of all post-masked bases")
        parser.add_argument('-o', '--output', metavar='OUTPUT_FILE',
                            default=sys.stdout,
                            help="Output file for the aligned FASTQ")
        args = parser.parse_args()
        self.__dict__.update( vars(args) )

    def initializeFromCall(self, fastqFile, outputFile, minQv):
        self.fastq = fastqFile
        logging.info('Creating a QualityMasker for "%s"' % self.fastq)
        # If no output file is set, default to STDOUT
        if outputFile is None:
            self.output = sys.stdout
        else:
            self.output = outputFile
        # If no minimum QV is set, use the default value
        if minQv is None:
            logging.info('No minimum QV specified, using default value %s' % self.minQV)
            self.minQv = MIN_QV
        else:
            logging.info('Minimum QV specified as %s' % minQv)
            self.minQv = minQv
        logging.info('No log-file set for this process')

    def validateSettings(self):
        filename, ext = os.path.splitext( self.fastq )
        try: # Try/Except Block for FASTQ input
            assert ext in ['.fq', '.fastq']
        except:
            raise ValueError("'%s' is not a recognized FASTQ file!" % self.fastq)
        #try:
        #    assert self.minQv > 0
        #except:
        #    raise ValueError("Minimum QV must be > 0!")

    ####################
    # Instance Methods #
    ####################

    def parseFastqData(self):
        self.fastqData = []
        logging.info('Reading Fastq data into memory from %s...' % self.fastq)
        for fastqRecord in FastqReader( self.fastq ):
            self.fastqData.append( fastqRecord )

    def trimFastqRecord(self, fastqRecord):
        # First identify any 5' N bases to trim
        if fastqRecord.sequence.startswith('N'):
            start = len( re.findall('^N+', fastqRecord.sequence)[0] )
        else:
            start = 0
        # Second identify any 3' N bases to trim
        if fastqRecord.sequence.endswith('N'):
            end = len(fastqRecord.sequence) - len( re.findall('N+$', fastqRecord.sequence)[0] )
        else:
            end = len(fastqRecord.sequence)
        # Third trim any non-N bases of low quality on the ends
        if self.minQv is not None:
            pass
        # Finally, subset the Fastq to the identified positions
        if start == 0 and end == len(fastqRecord.sequence):
            return fastqRecord
        else:
            return FastqRecord(fastqRecord.name,
                               fastqRecord.sequence[start:end],
                               qualityString=fastqRecord.qualityString[start:end])

    def trimFastqData(self):
        self.trimmedFastqs = []
        logging.info('Trimming low quality bases...')
        for fastqRecord in self.fastqData:
            trimmedFastq = self.trimFastqRecord( fastqRecord )
            self.trimmedFastqs.append( trimmedFastq )

    def filterFastqData(self):
        self.filteredFastqs = []
        logging.info('Filtering out short and empty sequences...')
        for fastqRecord in self.trimmedFastqs:
            if len(fastqRecord.sequence) < self.minLength:
                continue
            self.filteredFastqs.append( fastqRecord )

    def writeFastqData(self):
        logging.info('Writing the trimmed FASTQ data out to "%s"...' % self.output)
        with FastqWriter( self.output ) as writer: 
            for fastqRecord in self.filteredFastqs:
                writer.writeRecord( fastqRecord )

    def __call__(self):
        self.parseFastqData()
        self.trimFastqData()
        self.filterFastqData()
        self.writeFastqData()
        return self.output

if __name__ == '__main__':
    trimmer = QualityTrimmer()
    trimmer()
