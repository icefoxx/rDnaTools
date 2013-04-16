import sys
import os
import logging

from numpy import where

from ..io.FastqIO import FastqReader, FastqWriter

__version__ = "0.1"

class QualityMasker(object):
    """
    Tool for masking low-quality bases in FASTQ files 
    """
    MIN_QV = 15
    MIN_BASES = 1000

    ##########################
    # Initialization Methods #
    ##########################

    def __init__(self, fastqFile=None, outputFile=None, minQv=None):
        if fastqFile is None:
            self.initializeFromArgs()
        else:
            self.initializeFromCall(fastqFile, outputFile, minQv)
        self.validateSettings()

    def initializeFromArgs(self):
        args = self.parseArguments()
        self.__dict__.update( vars(args) )

    def parseArguments(self):
        import argparse
        desc = 'A tool for aligning FASTQ quality values to FASTA MSAs'
        parser = argparse.ArgumentParser(description=desc)
        parser.add_argument('fastq', metavar='FASTQ_FILE',
                            help="File of FASTQ sequence data to align")
        parser.add_argument('-q', '--minimum_qv', metavar='MIN_QV',
                            type=int, default=self.MIN_QV, dest='minQv',
                            help="Minimum QV to leave un-masked")
        parser.add_argument('-b', '--minimum_bases', metavar='MIN_BASES',
                            type=int, default=self.MIN_BASES, dest='minBases',
                            help="Minimum number of bases to require " + \
                                 "of all post-masked bases")
        parser.add_argument('-o', '--output', metavar='OUTPUT_FILE',
                            default=sys.stdout,
                            help="Output file for the aligned FASTQ")
        return parser.parse_args()

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
            self.minQv = self.MIN_QV
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
        try:
            assert self.minQv > 0
        except:
            raise ValueError("Minimum QV must be > 0!")

    ####################
    # Instance Methods #
    ####################

    def maskFastqRecord(self, fastqRecord):
        indicesToMask = where(fastqRecord.quality < self.minQv)[0]
        maskedSequence = ['-' if any(indicesToMask == pos) else base
                              for pos, base in enumerate(fastqRecord.sequence)]
        fastqRecord.sequence = ''.join(maskedSequence)
        return fastqRecord

    def run(self):
        self.parseFastqData()
        self.maskFastqData()
        self.writeFastqData()
        return self.output

    def parseFastqData(self):
        self.fastqData = []
        logging.info('Reading Fastq data into memory from %s...' % self.fastq)
        for fastqRecord in FastqReader( self.fastq ):
            self.fastqData.append( fastqRecord )

    def maskFastqData(self):
        self.maskedFastqs = []
        logging.info('Masking low quality bases')
        for fastqRecord in self.fastqData:
            maskedFastq = self.maskFastqRecord( fastqRecord )
            self.maskedFastqs.append( maskedFastq )

    def writeFastqData(self):
        logging.info('Writing the masked Fastq data out to "%s"...' % self.output)
        with FastqWriter( self.output ) as writer: 
            for fastqRecord in self.maskedFastqs:
                writer.writeRecord( fastqRecord )

if __name__ == '__main__':
    masker = QualityMasker()
    masker.run()
