#!/usr/bin/env python

import os
import sys
import logging
import subprocess

from pbtools.pbrdna._utils import which, fileExists, createDirectory
from pbtools.pbrdna.io import BasH5Extractor
from pbtools.pbrdna.io import SummaryReader
from pbtools.pbrdna.fastq import QualityAligner
from pbtools.pbrdna.fastq import QualityMasker
from pbtools.pbrdna.mothur import MothurFactory

__version__ = "0.1"

class rDnaPipeline(object):
    """
    A tool for running a community analysis pipeline on PacBioData
    """
    #####################
    # Variable Defaults #
    #####################

    MIN_DIST = 0.001
    MAX_DIST = 0.5
    MIN_QV = 15
    DEFAULT_FRAC = 0.8
    MIN_LENGTH = 500
    MIN_RATIO = 0.5
    CLUSTER_METHODS = ('nearest', 'average', 'furthest')

    ##########################
    # Initialization Methods #
    ##########################

    def __init__(self, sequenceFile=None):
        if sequenceFile is None:
            self.initializeFromArgs()
        else:
            self.initializeFromCall(sequenceFile)
        self.validateSettings()
        self.initializeOutput()
        self.initializeLogger()

    def initializeFromArgs(self):
        import argparse
        desc = 'A pipeline tool for analyzing rRNA amplicons'
        parser = argparse.ArgumentParser(description=desc)
        parser.add_argument('sequenceFile', metavar='FILE',
                            help="File of rRNA sequencing data to use")
        parser.add_argument('-d', '--distance', metavar='FLOAT', 
                            type=float, default=0.03,
                            help="Distance at which to cluster sequences")
        parser.add_argument('-n', '--num_processes', metavar='INT',
                            default=1, dest='numProc', type=int,
                            help="Number of processors to use")
        parser.add_argument('-f', '--fraction', metavar='FLOAT', 
                            type=float, default=self.DEFAULT_FRAC,
                            help='Fraction of full-length to require of each read')
        parser.add_argument('-c', '--clustering_method', metavar='METHOD',
                            dest='clusteringMethod', default='average',
                            choices=self.CLUSTER_METHODS,
                            help="Distance algorithm to use in clustering")
        parser.add_argument('-o', '--output', dest='outputDir', metavar='DIR',
                            default='rna_pipeline_run',
                            help="Specify the output folder")
        parser.add_argument('-q', '--minimum_qv', type=int, metavar='INT', 
                            dest='minQv', default=self.MIN_QV,
                            help='Minimun QV to allow after sequence masking')
        parser.add_argument('-l', '--minimum_length', type=int, metavar='INT', 
                            dest='minLength', default=self.MIN_LENGTH,
                            help='Minimun length sequence to allow after masking')
        parser.add_argument('-r', '--minimum_ratio', type=float, metavar='FLOAT',
                            dest='minRatio', default=self.MIN_RATIO,
                            help='Minimum ratio of retained bases to allow after masking')
        parser.add_argument('--disable_masking', action='store_true',
                            dest='disableMasking',
                            help="Turn off the low-quality masking step")
        parser.add_argument('--disable_resequencing', action='store_true',
                            dest='disableResequencing',
                            help="Turn off the resequencing step")
        parser.add_argument('-A', '--alignment_reference', metavar='REF',
                            default='silva.both.align', dest='alignmentRef',
                            help="Reference MSA for aligning query sequences")
        parser.add_argument('-C', '--chimera_reference', metavar='REF',
                            default='silva.gold.align', dest='chimeraRef',
                            help="Reference MSA for Chimera detection")
        parser.add_argument('--blasr', metavar='BLASR_PATH', 
                            help="Specify the path to the Blasr executable")
        parser.add_argument('--mothur', metavar='MOTHUR_PATH', 
                            help="Specify the path to the Mothur executable")
        parser.add_argument('--debug', action='store_true',
                            help="Turn on DEBUG message logging")
        args = parser.parse_args()
        self.__dict__.update( vars(args) )

    def validateSettings(self):
        # Searching for Mothur executable, and set the Mothur Process counter
        self.mothur = which('mothur')
        if self.mothur is None:
           raise OSError('Mothur executable not found!')
        self.processCount = 0
        # Validate the input file
        root, ext = self.splitRootFromExt( self.sequenceFile )
        if ext in ['.bas.h5', '.fofn']:
            self.dataType = 'bash5'
        elif ext in ['.fq', '.fastq']:
            self.dataType = 'fastq'
            self.disableResequencing = True
        elif ext in ['.fa', '.fsa', '.fasta']:
            self.dataType = 'fasta'
            self.disableMasking = True
            self.disableResequencing = True
        else:
            raise TypeError('Sequence file must be a bas.h5 file, a ' + \
                            'fasta file, or a fofn of multiple such files')
        # Validate the Num_Processes argument
        if self.numProc <= 0:
            raise ValueError("Number of processes must be >= 1!")
        # Validate the Distance argument
        try:
            dist = float(self.distance)
        except TypeError:
            raise TypeError("Distance is not a valid Float!")
        if dist < self.MIN_DIST: 
            raise ValueError("Distance must be > %s)" % self.MIN_DIST)
        if dist > self.MAX_DIST:
            raise ValueError("Distance must be < %s!" % self.MAX_DIST)

    def initializeOutput(self):
        # Create the Output directory
        createDirectory( self.outputDir )
        # Create a symbolic link from the data file to the output dir
        baseName = os.path.basename( self.sequenceFile )
        symlinkPath = os.path.join( self.outputDir, baseName )
        if os.path.exists( symlinkPath ):
            pass
        else:
            absPath = os.path.abspath( self.sequenceFile )
            os.symlink( absPath, symlinkPath )
        self.sequenceFile = baseName
        # Move into the Output directory and create Log directory and files
        os.chdir( self.outputDir )
        createDirectory( 'log' )
        stdoutLog = os.path.join('log', 'mothur_stdout.log')
        stderrLog = os.path.join('log', 'mothur_stderr.log')
        self.logFile = os.path.join('log', 'rna_pipeline.log')
        # Instantiate the MothurCommandFactory object
        self.factory = MothurFactory( self.mothur, 
                                      self.numProc, 
                                      stdoutLog, 
                                      stderrLog)

    def initializeLogger(self):
        dateFormat = "%Y-%m-%d %I:%M:%S"
        self.log = logging.getLogger()
        if self.debug:
            self.log.setLevel( logging.DEBUG )
        else:
            self.log.setLevel( logging.INFO )
        # Initialize the LogHandler for the master log file
        logHandler = logging.FileHandler( self.logFile )
        lineFormat = "%(asctime)s %(levelname)s %(processName)s " + \
                     "%(funcName)s %(lineno)d %(message)s"
        logFormatter = logging.Formatter(fmt=lineFormat, datefmt=dateFormat)
        logHandler.setFormatter( logFormatter )
        self.log.addHandler( logHandler )
        # Initialize a LogHandler for STDOUT
        outHandler = logging.StreamHandler( stream=sys.stdout )
        outLineFormat = "%(asctime)s %(message)s"
        outFormatter = logging.Formatter(fmt=outLineFormat, datefmt=dateFormat)
        outHandler.setFormatter( outFormatter )
        self.log.addHandler( outHandler )
        # Record the initialization of the pipeline
        self.log.info("INFO logger initialized")
        self.log.debug("DEBUG logger initialized")
        self.log.info("Initializing RnaPipeline v%s" % __version__)
        self.log.debug("Using the following parameters:")
        for param, value in self.__dict__.iteritems():
            self.log.debug("\t%s = %s" % (param, value))
        self.log.info("Initialization of RnaPipeline completed\n")

    #################
    # Class Methods #
    #################

    @classmethod
    def splitRootFromExt(cls, inputFile):
        root, ext = os.path.splitext( inputFile )
        if ext == '.h5':
            root, ext = os.path.splitext( root )
            return (root, ext+'.h5')
        return (root, ext)

    @classmethod
    def predictOutputFile(cls, inputFile, outputType):
        root, ext = cls.splitRootFromExt( inputFile )
        return root + '.' + outputType

    ####################
    # Instance Methods #
    ####################

    def getProcessLogFile(self, process, isMothurProcess=False):
        if isMothurProcess:
            logFile = 'process%02d.mothur.%s.logfile' % (self.processCount, 
                                                         process)
        else:
            logFile = 'process%02d.%s.logfile' % (self.processCount, process)
        return os.path.join('log', logFile)

    def processSetup(self, inputFile, processName, outputSuffix):
        """ 
        Return a tuple containing the output file and a boolean flag describing
        whether the output file already exists
        """
        self.log.info('Preparing to run %s on "%s"' % (processName, inputFile))
        self.processCount += 1
        outputFile = self.predictOutputFile(inputFile, outputSuffix)
        if fileExists( outputFile ):
            self.log.info('Output file "%s" detected, skipping process...\n' % \
                                                                     outputFile)
        else:
            self.log.info('Output file "%s" not detected, running process...' % \
                                                                      outputFile)
        return outputFile

    def processCleanup(self, outputFile):
        """
        Log if the process successfully created it's output, and raise an
        error message if not
        """
        if fileExists( outputFile ):
            self.log.info('"%s" created - process completed successfully\n' % \
                                                                    outputFile)
        else:
            self.log.info('"%s" not created - process failure!\n' % outputFile)
            raise IOError('"%s" not created - process failure!' % outputFile)

    def convertFastqToFasta(self, fastqFile):
        outputFile = self.processSetup(fastqFile, 'Fastq.Info', 'fasta')
        if not fileExists( outputFile ):
            mothurArgs = {'fastq':fastqFile, 'fasta':'T', 'qfile':'F'}
            logFile = self.getProcessLogFile('fastq.info', True)
            self.factory.runJob('fastq.info', mothurArgs, logFile)
            self.processCleanup( outputFile )
        return outputFile

    def alignSequences(self, fastaFile):
        outputFile = self.processSetup(fastaFile, 'Align.Seqs', 'align')
        if not fileExists( outputFile ):
            mothurArgs = {'fasta':fastaFile,
                          'reference':self.alignmentRef,
                          'flip':'t'}
            logFile = self.getProcessLogFile('align.seqs', True)
            self.factory.runJob('align.seqs', mothurArgs, logFile)
            self.processCleanup( outputFile )
        return outputFile

    def screenSequences(self, alignFile, start=None, end=None, minLength=None):
        if alignFile.endswith('.align'):
            outputExt = 'good.align'
        elif alignFile.endswith('.fasta'):
            outputExt = 'good.fasta'
        outputFile = self.processSetup(alignFile, 'Screen.Seqs', outputExt)
        if not fileExists( outputFile ):
            mothurArgs = {'fasta':alignFile,
                          'start':start,
                          'end':end,
                          'minlength':minLength}
            logFile = self.getProcessLogFile('screen.seqs', True)
            self.factory.runJob('screen.seqs', mothurArgs, logFile)
            self.processCleanup( outputFile )
        return outputFile

    def summarizeSequences(self, fastaFile):
        outputFile = self.processSetup(fastaFile, 'Summary.Seqs', 'summary')
        if not fileExists( outputFile ):
            mothurArgs = {'fasta':fastaFile}
            logFile = self.getProcessLogFile('summary.seqs', True)
            self.factory.runJob('summary.seqs', mothurArgs, logFile)
            self.processCleanup( outputFile )
        return outputFile

    def parseSummaryFile(self, summaryFile):
        self.log.info('Preparing to run SummaryParser...')
        parser = SummaryParser(summaryFile, self.fraction)
        self.log.info('Identifying full-length alignment positions...')
        start, end = parser.getFullLengthPositions()
        self.log.info('Full-length start is NAST Alignment position %s' % start)
        self.log.info('Full-length end is NAST Alignment position %s' % end)
        self.log.info('Calculating minimum allowed alignment positions...')
        maxStart, minEnd = parser.getAllowedPositions()
        self.log.info('Maximum allowed start is NAST Alignment position %s' % maxStart)
        self.log.info('Minimum allowed end is NAST Alignment position %s\n' % minEnd)
        return maxStart, minEnd

    def findChimeras(self, alignFile):
        outputFile = self.processSetup(alignFile, 'UCHIME', 'uchime.accnos')
        if not fileExists( outputFile ):
            mothurArgs = {'fasta':alignFile,
                          'reference':self.chimeraRef}
            logFile = self.getProcessLogFile('chimera.uchime', True)
            self.factory.runJob('chimera.uchime', mothurArgs, logFile)
            self.processCleanup( outputFile )
        return outputFile

    def removeSequences(self, alignFile, idFile):
        outputFile = self.processSetup(alignFile, 'Remove.Seqs', 'pick.align')
        if not fileExists( outputFile ):
            mothurArgs = {'fasta':alignFile,
                          'accnos':idFile}
            logFile = self.getProcessLogFile('remove.seqs', True)
            self.factory.runJob('remove.seqs', mothurArgs, logFile)
            self.processCleanup( outputFile )
        return outputFile

    def filterSequences(self, alignFile):
        outputFile = self.processSetup(alignFile, 'Filter.Seqs', 'filter.fasta')
        if not fileExists( outputFile ):
            mothurArgs = {'fasta':alignFile,
                          'vertical':'T'}
            logFile = self.getProcessLogFile('filter.seqs', True)
            self.factory.runJob('filter.seqs', mothurArgs, logFile)
            self.processCleanup(outputFile)
        return outputFile

    def extractCcsFromBasH5(self, inputFile):
        outputFile = self.processSetup(inputFile, 'extractCcsFromBasH5', 'fastq')
        if not fileExists( outputFile ):
            extractor = BasH5Extractor(inputFile, outputFile)
            extractor.outputCcsFastq()
            self.processCleanup(outputFile)
        return outputFile

    def addQualityToAlignment(self, fastqFile, alignFile):
        outputFile = self.processSetup(alignFile, 'FastqAligner', 'fastq')
        if not fileExists( outputFile ):
            aligner = FastqAligner(fastqFile, alignFile, outputFile)
            aligner.run()
            self.processCleanup(outputFile)
        return outputFile

    def maskFastqSequences(self, fastqFile):
        outputFile = self.processSetup(fastqFile, 'FastqMasker', 'masked.fastq')
        if not fileExists( outputFile ):
            masker = FastqMasker(fastqFile, outputFile, self.minQv)
            masker.run()
            self.processCleanup(outputFile)
        return outputFile

    def calculateDistanceMatrix(self, alignFile):
        outputFile = self.processSetup(alignFile, 'Dist.Seqs', 'phylip.dist')
        if not fileExists( outputFile ):
            mothurArgs = {'fasta':alignFile,
                          'calc':'nogaps',
                          'output':'lt'}
            logFile = self.getProcessLogFile('dist.seqs', True)
            self.factory.runJob('dist.seqs', mothurArgs, logFile)
            self.processCleanup(outputFile)
        return outputFile

    def clusterSequences(self, distanceMatrix):
        if self.clusteringMethod == 'nearest':
            outputSuffix = 'nn.list'
        elif self.clusteringMethod == 'average':
            outputSuffix = 'an.list'
        elif self.clusteringMethod == 'furthest':
            outputSuffix = 'fn.list'
        outputFile = self.processSetup(distanceMatrix, 'Cluster', outputSuffix)
        if not fileExists( outputFile ):
            mothurArgs = {'phylip':distanceMatrix,
                          'method':self.clusteringMethod}
            logFile = self.getProcessLogFile('cluster', True)
            self.factory.runJob('cluster', mothurArgs, logFile)
            self.processCleanup(outputFile)
        return outputFile

    def __call__(self):
        # Extract and convert the data to FASTA format as needed
        if self.dataType == 'bash5':
            basH5File = self.sequenceFile
            fastqFile = self.extractCcsFromBasH5( basH5File )
            fastaFile = self.convertFastqToFasta( fastqFile )
        elif self.dataType == 'fastq':
            fastqFile = self.sequenceFile
            fastaFile = self.convertFastqToFasta( fastqFile )
        elif self.dataType == 'fasta':
            fastaFile = self.sequenceFile
        # Align the FASTA sequences via Mothur's NAST implementation
        alignedFile = self.alignSequences( fastaFile )
        # Summarize the start and end positions, and sequence lengths
        summaryFile = self.summarizeSequences( alignedFile )
        maxStart, minEnd = self.parseSummaryFile( summaryFile )
        # Remove sequences that don't cover the appropriate regions 
        screenedFile = self.screenSequences(alignedFile, 
                                            start=maxStart,
                                            end=minEnd)
        # Identify and remove chimeric reads
        chimeraIds = self.findChimeras( screenedFile )
        noChimeraFile = self.removeSequences( screenedFile, chimeraIds )
        # Filter out un-used columns to speed up re-alignment and clustering
        filteredFile = self.filterSequences( noChimeraFile )
        # If Masking is disabled, use the filtered sequence file for
        # Clustering and skip to that step
        if self.disableMasking:
            fileForClustering = filteredFile
        # Otherwise, create an aligned FASTQ, mask the low-quality bases,
        #   and emove any sequences that have to few bases remaining
        if not self.disableMasking:
            alignedFastqFile = self.addQualityToAlignment( fastqFile, filteredFile )
            maskedFastq = self.maskFastqSequences( alignedFastqFile )
            maskedFasta = self.convertFastqToFasta( maskedFastq )
            screenedFasta = self.screenSequences(maskedFasta,
                                                 minLength=self.minLength)
            fileForClustering = screenedFasta
        # Calculate sequence distances and cluster
        distanceMatrix = self.calculateDistanceMatrix( fileForClustering )
        self.clusterSequences( distanceMatrix )


if __name__ == '__main__':
    rdnap = rDnaPipeline()
    rdnap()
