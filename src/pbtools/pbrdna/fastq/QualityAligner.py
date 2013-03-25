#!/usr/bin/env python

import sys
import re
import os
import csv
import subprocess
import logging

from collections import namedtuple
from random import random
from string import maketrans

from ..io.FastaIO import FastaRecord, FastaReader, FastaWriter
from ..io.FastqIO import FastqRecord, FastqReader, FastqWriter

__version__ = "0.1"

class FastqAligner(object):
    """
    Tool for aligning quality values from FASTQ files to trimmed and  
    gapped multi-sequence alignments of their sequences
    """

    DNA_TRANSLATOR = maketrans('AGCT', 'TCGA')

    BlasrRecord = namedtuple('BlasrRecord', ['qname', 'tname', 'qstrand', 'tstrand',
                                             'score', 'pctsimilarity', 
                                             'tstart', 'tend', 'tlength',
                                             'qtsart', 'qend', 'qlength', 'ncells'] )
    ##########################
    # Initialization Methods #
    ##########################

    def __init__(self, fastqFile=None, alignedFile=None, outputFile=None):
        if fastqFile is None or alignedFile is None:
            self.initializeFromArgs()
        else:
            self.initializeFromCall(fastqFile, alignedFile, outputFile)
        self.validateSettings()

    def initializeFromArgs(self):
        args = self.parseArguments()
        self.__dict__.update( vars(args) )

    def parseArguments(self):
        import argparse
        desc = 'A tool for  aligning FASTQ quality values to FASTA MSAs'
        parser = argparse.ArgumentParser(description=desc)
        parser.add_argument('fastq', metavar='FASTQ_FILE',
                            help="File of FASTQ sequence data to align")
        parser.add_argument('aligned', metavar='ALIGNMENT_FILE',
                            help="Gapped FASTA sequences to align to")
        parser.add_argument('-b', '--blasr', metavar='BLASR_PATH',
                            help="Path to the local Blasr executable")
        parser.add_argument('-o', '--output', metavar='OUTPUT_FILE',
                            default=sys.stdout,
                            help="Output file for the aligned FASTQ")
        return parser.parse_args()

    def initializeFromCall(self, fastqFile, alignedFile, outputFile):
        self.fastq = fastqFile
        self.aligned = alignedFile
        if outputFile is None:
            self.output = sys.stdout
        else:
            self.output = outputFile

    def validateSettings(self):
        filename, ext = os.path.splitext( self.fastq )
        logging.info('Pulling QV data from "%s"' % self.fastq)
        try: # Try/Except Block for FASTQ input
            assert ext in ['.fq', '.fastq']
        except:
            raise ValueError("'%s' is not a recognized FASTQ file!" % self.fastq)
        filename, ext = os.path.splitext( self.aligned )
        try: # Try/Except Block for Alignment input
            assert ext in ['.fa', '.fsa', '.fasta', '.align']
        except:
            raise ValueError("'%s' is not a recognized FASTA file!" % self.aligned)
        logging.info('Creating a FastqAligner for "%s"' % self.aligned)
        logging.info('No log-file set for this process')

    #################
    # Class Methods #
    #################

    @classmethod
    def which(cls, program):
        """
        Find and return path to local executables  
        """
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file
        return None

    @classmethod
    def convertFastqToFasta(cls, fastqRecord):
        return FastaRecord( fastqRecord.name, fastqRecord.sequence )

    @classmethod
    def reverseComplement(cls, record):
        if isinstance(record, str):
            return record[::-1].translate( cls.DNA_TRANSLATOR )
        elif isinstance(record, FastaRecord):
            return FastaRecord(record.name, 
                               cls.reverseComplement( record.sequence ))
        elif isinstance(record, FastqRecord):
            return FastqRecord(record.name, 
                               cls.reverseComplement( record.sequence ),
                               record.quality[::-1])
        else:
            raise ValueError("Record must be either FASTA or FASTQ")

    @classmethod
    def getZmw(cls, record):
        parts = record.name.split('/')
        zmw = '/'.join( parts[0:2] )
        return zmw

    @classmethod
    def getSeqParts(cls, record):
        parts = re.split('(\.+|-+|\w+)', record.sequence)
        parts = [len(p) if re.match('(\.|-)', p) else p
                        for p in parts 
                        if p != '']
        return parts

    @classmethod
    def createUnalignedRecord(cls, seqParts, zmw):
        bases = [b for b in seqParts if type(b) is str]
        unalignedSequence = ''.join(bases)
        unalignedRecord = FastaRecord( zmw, unalignedSequence )
        return unalignedRecord

    @classmethod
    def runBlasr(cls, fastqRecord, alignedRecord):
        # Write the query and reference records to file
        tempId = str(int(random() * 10000000))
        tempRef = 'temp_ref_%s.fasta' % tempId
        with FastaWriter(tempRef) as handle:
            fastaRecord = cls.convertFastqToFasta( fastqRecord )
            handle.writeRecord( fastaRecord )
        tempQuery = 'temp_query_%s.fasta' % tempId 
        with FastaWriter(tempQuery) as handle:
            handle.writeRecord( alignedRecord )
        # Create and run the command-line
        tempOut = 'temp_%s.m1' % tempId
        cline = 'blasr %s %s -m 1 -bestn 1 -out %s' % (tempQuery,
                                                       tempRef,
                                                       tempOut)
        p = subprocess.Popen( cline.split() )
        stdout, stderr = p.communicate()
        # Parse and return the best hit and remove temp files
        bestHit = cls.readBestBlasrHit(tempOut)
        os.remove(tempRef)
        os.remove(tempQuery)
        os.remove(tempOut)
        return bestHit

    @classmethod
    def readBestBlasrHit(cls, blasrOutput):
        with open(blasrOutput, 'r') as handle:
            csvReader = csv.reader(handle, delimiter=' ')
            blasrHits = list( map(cls.BlasrRecord._make, csvReader) )
        return blasrHits[0]

    @classmethod
    def trimFastqRecord(cls, fastqRecord, blasrHit):
        if blasrHit.qstrand != blasrHit.tstrand:
            fastqRecord = cls.reverseComplement( fastqRecord )
        start = int(blasrHit.tstart)
        end = int(blasrHit.tend)
        trimmedSequence = fastqRecord.sequence[start:end]
        trimmedQualities = fastqRecord.quality[start:end]
        trimmedRecord = FastqRecord(fastqRecord.name,
                                    trimmedSequence,
                                    trimmedQualities)
        return trimmedRecord

    @classmethod
    def addGappedQualities(cls, fastqRecord, seqParts, alignedRecord):
        qualities = ''
        pos = 0
        for part in seqParts:
            if type(part) is int:
                qualities += '!' * part
            else:
                qualities += fastqRecord.qualityString[pos:pos+len(part)]
                pos += len(part)
        alignedRecord = FastqRecord(alignedRecord.name,
                                    alignedRecord.sequence,
                                    qualityString=qualities) 
        return alignedRecord

    ####################
    # Instance Methods #
    ####################

    def run(self):
        self.parseFastqData()
        self.alignFastqData()
        self.writeFastqData()
        return self.output

    def parseFastqData(self):
        self.sequenceData = {}
        logging.info('Reading QV data from "%s"...' % self.fastq)
        counter = 0
        for record in FastqReader( self.fastq ):
            zmw = self.getZmw( record )
            self.sequenceData[zmw] = record
            counter += 1
        logging.info('A total of %s Fastq records were read into memory' % counter)

    def alignFastqData(self):
        self.alignedFastqs = []
        logging.info('Combining data from the Aligned Fasta and Fastq files...')
        counter = 0
        for record in FastaReader( self.aligned ):
            zmw = self.getZmw( record )
            try:
                fastqRecord = self.sequenceData[zmw]
            except KeyError: 
                raise KeyError("No quality data found for '%s'!" % zmw)
            seqParts = self.getSeqParts(record)
            unalignedRecord = self.createUnalignedRecord(seqParts, zmw)
            blasrHit = self.runBlasr(fastqRecord, unalignedRecord)
            trimmedFastq = self.trimFastqRecord(fastqRecord, blasrHit)
            try:
                assert unalignedRecord.sequence == trimmedFastq.sequence
            except AssertionError:
                raise ValueError("Sequences don't match for '%s'" % zmw)
            updatedRecord = self.addGappedQualities(trimmedFastq, seqParts, record)
            counter += 1
            self.alignedFastqs.append(updatedRecord)
        logging.info('A total of %s aligned Fastq records were created' % counter)

    def writeFastqData(self):
        logging.info('Writing aligned Fastq data out to "%s"' % self.output)
        with FastqWriter( self.output ) as handle:
            for alignedFastq in self.alignedFastqs:
                handle.writeRecord( alignedFastq )


if __name__ == '__main__':
    aligner = FastqAligner()
    aligner.run()
