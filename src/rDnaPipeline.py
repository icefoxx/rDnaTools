#!/home/UNIXHOME/bbowman/environments/rDnaTools/bin/python2.7

import os
import sys
import logging
import subprocess

from pbrdna._utils import *
from pbrdna.io.BasH5IO import BasH5Extractor
from pbrdna.io.MothurIO import SummaryReader
from pbrdna.fastq import QualityFilter
from pbrdna.fastq.QualityAligner import QualityAligner
from pbrdna.fastq.QualityMasker import QualityMasker
from pbrdna.mothur.MothurTools import MothurRunner
from pbrdna.cluster.ClusterSeparator import ClusterSeparator
from pbrdna.resequence.DagConTools import DagConRunner

__version__ = "0.2"

MIN_DIST = 0.001
MAX_DIST = 0.5
MIN_ACCURACY = 0.99
MIN_QV = 15
DEFAULT_FRAC = 0.8
MIN_LENGTH = 500
MIN_RATIO = 0.5
PRECLUSTER_DIFFS = 4
CLUSTER_METHODS = ('nearest', 'average', 'furthest')
DEFAULT_METHOD = 'average'

class rDnaPipeline( object ):
    """
    A tool for running a community analysis pipeline on PacBioData
    """

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
        parser.add_argument('-a', '--minimum_accuracy', type=float, metavar='FLOAT',
                            dest='minAccuracy', default=MIN_ACCURACY,
                            help='Minimum predicted sequence accuracy')
        parser.add_argument('-d', '--distance', metavar='FLOAT', 
                            type=float, default=0.03,
                            help="Distance at which to cluster sequences")
        parser.add_argument('-n', '--num_processes', metavar='INT',
                            default=1, dest='numProc', type=int,
                            help="Number of processors to use")
        parser.add_argument('-f', '--fraction', metavar='FLOAT', 
                            type=float, default=DEFAULT_FRAC,
                            help='Fraction of full-length to require of each read')
        parser.add_argument('-c', '--clustering_method', metavar='METHOD',
                            dest='clusteringMethod', default=DEFAULT_METHOD,
                            choices=CLUSTER_METHODS,
                            help="Distance algorithm to use in clustering")
        parser.add_argument('-o', '--output', dest='outputDir', metavar='DIR',
                            default='rna_pipeline_run',
                            help="Specify the output folder")
        parser.add_argument('-q', '--minimum_qv', type=int, metavar='INT', 
                            dest='minQv', default=MIN_QV,
                            help='Minimum QV to allow after sequence masking')
        parser.add_argument('-l', '--minimum_length', type=int, metavar='INT', 
                            dest='minLength', default=MIN_LENGTH,
                            help='Minimun length sequence to allow after masking')
        parser.add_argument('--precluster_diffs', type=int, metavar='INT',
                            dest='preclusterDiffs', default=PRECLUSTER_DIFFS,
                            help='Maximum number of differences to allow in pre-clustering')
        parser.add_argument('-r', '--minimum_ratio', type=float, metavar='FLOAT',
                            dest='minRatio', default=MIN_RATIO,
                            help='Minimum ratio of retained bases to allow after masking')
        parser.add_argument('-A', '--alignment_reference', metavar='REF',
                            default='silva.both.align', dest='alignmentRef',
                            help="Reference MSA for aligning query sequences")
        parser.add_argument('-C', '--chimera_reference', metavar='REF',
                            default='silva.gold.align', dest='chimeraRef',
                            help="Reference MSA for Chimera detection")
        parser.add_argument('--enable_masking', action='store_true',
                            dest='enableMasking',
                            help="Turn off the low-quality Masking step")
        parser.add_argument('--disable_clustering', action='store_false',
                            dest='enableClustering',
                            help="Turn off the Clustering and Resequencing steps")
        parser.add_argument('--disable_consensus', action='store_false',
                            dest='enableConsensus',
                            help="Turn off the Consensus step")
        parser.add_argument('--blasr', metavar='BLASR_PATH', 
                            help="Specify the path to the Blasr executable")
        parser.add_argument('--mothur', metavar='MOTHUR_PATH', default='mothur',
                            help="Specify the path to the Mothur executable")
        parser.add_argument('--debug', action='store_true',
                            help="Turn on DEBUG message logging")
        args = parser.parse_args()
        self.__dict__.update( vars(args) )

    def validateSettings(self):
        # Validate the input file
        root, ext = self.splitRootFromExt( self.sequenceFile )
        if ext in ['.bas.h5', '.fofn']:
            self.dataType = 'bash5'
        elif ext in ['.fq', '.fastq']:
            self.dataType = 'fastq'
        elif ext in ['.fa', '.fsa', '.fasta']:
            self.dataType = 'fasta'
            self.enableMasking = False
            self.enableConsensus = False
        else:
            raise TypeError('Sequence file must be a bas.h5 file, a ' + \
                            'fasta file, or a fofn of multiple such files')
        # If Clustering was disabled, also disable the consensus process
        if not self.enableClustering:
            self.enableConsensus = False
        # If Consensus is enabled, initialize the appropriate tool
        if self.enableConsensus:
            self.consensusTool = DagConRunner('gcon.py', 'r')
        # Searching for Mothur executable, and set the Mothur Process counter
        self.mothur = validateExecutable( self.mothur )
        self.processCount = 0
        # Validate numerical parameters
        validateInt( self.numProc, minValue=0 )
        validateFloat( self.distance, minValue=MIN_DIST, maxValue=MAX_DIST )

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
        # Instantiate the MothurRunner object
        self.factory = MothurRunner( self.mothur, 
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

    def getProcessLogFile(self, process, isMothurProcess=False):
        if isMothurProcess:
            logFile = 'process%02d.mothur.%s.logfile' % (self.processCount, 
                                                         process)
        else:
            logFile = 'process%02d.%s.logfile' % (self.processCount, process)
        return os.path.join('log', logFile)

    def processSetup(self, inputFile, processName, suffix=None, suffixList=None):
        """ 
        Return a tuple containing the output file and a boolean flag describing
        whether the output file already exists
        """
        self.log.info('Preparing to run %s on "%s"' % (processName, inputFile))
        self.processCount += 1
        if suffix:
            outputFile = self.predictOutputFile(inputFile, suffix)
            return outputFile
        elif suffixList:
            outputFiles = []
            for suffix in suffixList:
                outputFile = self.predictOutputFile( inputFile, suffix )
                outputFiles.append( outputFile )
            return outputFiles

    def outputFilesExist( self, outputFile=None, outputList=None ):
        if outputFile:
            if fileExists( outputFile ):
                self.log.info('Output files detected, skipping process...\n')
                return True
            else:
                self.log.info('Output files not found, running process...')
                return False
        elif outputList:
            if allFilesExist( outputList ):
                self.log.info('Output files detected, skipping process...\n')
                return True
            else:
                self.log.info('Output files not found, running process...')
                return False

    def checkOutputFile( self, outputFile ):
        if fileExists( outputFile ):
            self.log.info('Expected output "%s" found' % outputFile)
        else:
            msg = 'Expected output "%s" not found!' % outputFile
            self.log.info( msg )
            raise IOError( msg )

    def processCleanup(self, outputFile=None, outputList=None):
        """
        Log if the process successfully created it's output, and raise an
        error message if not
        """
        if outputFile:
            self.checkOutputFile( outputFile )
        elif outputList:
            for outputFile in outputList:
                self.checkOutputFile( outputFile )
        self.log.info('All expected output files found - process successful!\n')

    def writeDummyFile(self, dummyFile):
        with open(dummyFile, 'w') as handle:
            handle.write('DONE')
        return dummyFile

    def extractCcsFromBasH5(self, inputFile):
        outputFile = self.processSetup( inputFile, 
                                        'extractCcsFromBasH5', 
                                        suffix='fastq' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        extractor = BasH5Extractor( inputFile, outputFile )
        extractor.outputCcsFastq()
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def filterFastqFile(self, fastqFile):
        outputList = self.processSetup( fastqFile, 
                                        'FilterQuality', 
                                        suffix='filter.fastq' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        aligner = QualityFilter( fastqFile, outputFile, self. )
        aligner()
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def separateFastqFile(self, fastqFile):
        outputList = self.processSetup( fastqFile, 
                                        'Fastq.Info', 
                                        suffixList=['fasta', 'qual'] )
        if self.outputFilesExist( outputList=outputList ):
            return outputList
        mothurArgs = {'fastq':fastqFile, 'fasta':'T', 'qfile':'T'}
        logFile = self.getProcessLogFile('fastq.info', True)
        self.factory.runJob('fastq.info', mothurArgs, logFile)
        self.processCleanup( outputList=outputList )
        return outputList

    def alignSequences(self, fastaFile):
        outputFile = self.processSetup( fastaFile, 
                                        'Align.Seqs', 
                                        suffix='align' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':fastaFile,
                      'reference':self.alignmentRef,
                      'flip':'t'}
        logFile = self.getProcessLogFile('align.seqs', True)
        self.factory.runJob('align.seqs', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def screenSequences(self, alignFile, start=None, end=None, minLength=None):
        if alignFile.endswith('.align'):
            outputExt = 'good.align'
        elif alignFile.endswith('.fasta'):
            outputExt = 'good.fasta'
        outputFile = self.processSetup( alignFile, 
                                         'Screen.Seqs', 
                                         suffix=outputExt )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'start':start,
                      'end':end,
                      'minlength':minLength}
        logFile = self.getProcessLogFile('screen.seqs', True)
        self.factory.runJob('screen.seqs', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def summarizeSequences(self, fastaFile):
        outputFile = self.processSetup( fastaFile, 
                                        'Summary.Seqs', 
                                        suffix='summary' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':fastaFile}
        logFile = self.getProcessLogFile('summary.seqs', True)
        self.factory.runJob('summary.seqs', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def parseSummaryFile(self, summaryFile):
        self.log.info('Preparing to run SummaryReader...')
        parser = SummaryReader(summaryFile, self.fraction)
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
        outputFile = self.processSetup( alignFile, 
                                        'UCHIME', 
                                        suffix='uchime.accnos' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'reference':self.chimeraRef}
        logFile = self.getProcessLogFile('chimera.uchime', True)
        self.factory.runJob('chimera.uchime', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def removeSequences(self, alignFile, idFile):
        outputFile = self.processSetup( alignFile, 
                                        'Remove.Seqs', 
                                        suffix='pick.align' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'accnos':idFile}
        logFile = self.getProcessLogFile('remove.seqs', True)
        self.factory.runJob('remove.seqs', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile
  
    def filterSequences(self, alignFile):
        outputFile = self.processSetup( alignFile, 
                                        'Filter.Seqs', 
                                        suffix='filter.fasta' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'fasta':alignFile,
                      'vertical':'T'}
        logFile = self.getProcessLogFile( 'filter.seqs', True )
        self.factory.runJob( 'filter.seqs', mothurArgs, logFile )
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def addQualityToAlignment(self, fastqFile, alignFile):
        outputFile = self.processSetup( alignFile, 
                                        'QualityAligner', 
                                        suffix='fastq' )
        if self.outputFilesExist( outputFile=output ):
            return output
        aligner = QualityAligner( fastqFile, alignFile, outputFile )
        aligner.run()
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def maskFastqSequences(self, fastqFile):
        outputFile = self.processSetup( fastqFile, 
                                        'QualityMasker', 
                                        suffix='masked.fastq' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        masker = QualityMasker(fastqFile, outputFile, self.minQv)
        masker.run()
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def uniqueSequences( self, alignFile ):
        if alignFile.endswith('.align'):
            outputSuffixes = ['unique.align', 'names']
        elif alignFile.endswith('.fasta'):
            outputSuffixes = ['unique.fasta', 'names']
        outputList = self.processSetup( alignFile,
                                        'Unique.Seqs',
                                        suffixList=outputSuffixes )
        if self.outputFilesExist( outputList=outputList ):
            return outputList
        mothurArgs = {'fasta':alignFile}
        logFile = self.getProcessLogFile('unique.seqs', True)
        self.factory.runJob('unique.seqs', mothurArgs, logFile)
        self.processCleanup( outputList=outputList )
        return outputList

    def preclusterSequences( self, alignFile, nameFile ):
        if alignFile.endswith('.align'):
            outputSuffixes = ['precluster.align', 'precluster.names']
        elif alignFile.endswith('.fasta'):
            outputSuffixes = ['precluster.fasta', 'precluster.names']
        outputList = self.processSetup( alignFile,
                                        'Pre.Cluster',
                                        suffixList=outputSuffixes )
        if self.outputFilesExist( outputList=outputList ):
            return outputList
        mothurArgs = { 'fasta':alignFile,
                       'name': nameFile,
                       'diffs':self.preclusterDiffs }
        logFile = self.getProcessLogFile('pre.cluster', True)
        self.factory.runJob('pre.cluster', mothurArgs, logFile)
        self.processCleanup( outputList=outputList )
        return outputList

    def calculateDistanceMatrix( self, alignFile ):
        outputFile = self.processSetup( alignFile, 
                                        'Dist.Seqs', 
                                        suffix='phylip.dist' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = { 'fasta':alignFile,
                       'calc':'nogaps',
                       'output':'lt' }
        logFile = self.getProcessLogFile('dist.seqs', True)
        self.factory.runJob('dist.seqs', mothurArgs, logFile)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def clusterSequences(self, distanceMatrix):
        if self.clusteringMethod == 'nearest':
            outputSuffix = 'nn.list'
        elif self.clusteringMethod == 'average':
            outputSuffix = 'an.list'
        elif self.clusteringMethod == 'furthest':
            outputSuffix = 'fn.list'
        outputFile = self.processSetup( distanceMatrix, 
                                        'Cluster', 
                                        suffix=outputSuffix )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        mothurArgs = {'phylip':distanceMatrix,
                      'method':self.clusteringMethod}
        logFile = self.getProcessLogFile( 'cluster', True )
        self.factory.runJob( 'cluster', mothurArgs, logFile )
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def separateClusterSequences(self, listFile, sequenceFile):
        outputFile = self.processSetup( listFile, 
                                        'ClusterSeparator', 
                                        suffix='list.clusters')
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        separator = ClusterSeparator( listFile, 
                                      sequenceFile, 
                                      self.distance, 
                                      outputFile )
        separator()
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def generateConsensusSequences(self, clusterListFile):
        outputFile = self.processSetup( clusterListFile, 
                                        'ClusterResequencer', 
                                        suffix='consensus')
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        consensusFiles = []
        with open( clusterListFile ) as handle:
            for line in handle:
                sequenceFile, referenceFile = line.strip().split()
                if referenceFile.endswith('None'):
                    consensusFiles.append( (sequenceFile, 'None') )
                else:
                    consensus = self.consensusTool( sequenceFile, referenceFile )
                    consensusFiles.append( (referenceFile, consensus) )
        with open( outputFile, 'w' ) as handle:
            for filenamePair in consensusFiles:
                handle.write('%s\t%s\n' % filenamePair)
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def cleanupConsensusFolder( self, consensusFile ):
        outputFile = self.processSetup( consensusFile, 
                                        'ConsensusCleanup', 
                                        suffix='consensus.cleanup' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        reseqPath = os.path.join( os.getcwd(), 'reseq' )
        for filename in os.listdir( reseqPath ):
            filePath = os.path.join( reseqPath, filename )
            if filePath.endswith('_input.fa'):
                os.remove( filePath )
            elif filePath.endswith('_input.fa.aln'):
                os.remove( filePath )
            elif filePath.endswith('_input.fa.aln_unsorted'):
                os.remove( filePath )
        self.writeDummyFile( outputFile )
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def selectFinalSequences( self, consensusFile ):
        outputFile = self.processSetup( consensusFile, 
                                        'SequenceSelector', 
                                        suffix='consensus.selected' )
        if self.outputFilesExist( outputFile=outputFile ):
            return outputFile
        selectedFiles = []
        with open( consensusFile ) as handle:
            for line in handle:
                referenceFile, consensusFile = line.strip().split()
                if consensusFile.endswith('None'):
                    selectedFiles.append( referenceFile )
                elif fasta_count( consensusFile ) == 1:
                    selectedFiles.append( consensusFile )
                else:
                    selectedFiles.append( referenceFile )
        with open( outputFile, 'w' ) as handle:
            for filename in selectedFiles:
                handle.write(filename + '\n')
        self.processCleanup( outputFile=outputFile )
        return outputFile

    def __call__(self):
        if self.dataType == 'bash5':
            fastqFile = self.extractCcsFromBasH5( self.sequenceFile )
        elif self.dataType == 'fastq':
            fastqFile = self.sequenceFile
        elif self.dataType == 'fasta':
            fastqFile = None
            fastaFile = self.sequenceFile
        # If we have a Fastq, filter low-quality reads and convert to FASTA
        if fastqFile:
            filteredFastq = self.filterFastqFile( fastqFile )
            fastaFile, qualFile = self.separateFastqFile( fastqFile )
"""
        # Align the Fasta sequences and remove partial reads
        alignedFile = self.alignSequences( fastaFile )
        summaryFile = self.summarizeSequences( alignedFile )
        maxStart, minEnd = self.parseSummaryFile( summaryFile )
        screenedFile = self.screenSequences(alignedFile, 
                                            start=maxStart,
                                            end=minEnd)
        # Identify and remove chimeric reads
        chimeraIds = self.findChimeras( screenedFile )
        noChimeraFile = self.removeSequences( screenedFile, chimeraIds )
        # Filter out un-used columns to speed up re-alignment and clustering
        filteredFile = self.filterSequences( noChimeraFile )
        # If masking is enabled, create an aligned FASTQ, mask the 
        # low-quality bases and remove over-masked reads
        if self.enableMasking:
            alignedFastqFile = self.addQualityToAlignment( fastqFile, filteredFile )
            maskedFastq = self.maskFastqSequences( alignedFastqFile )
            maskedFasta = self.convertFastqToFasta( maskedFastq )
            screenedFasta = self.screenSequences( maskedFasta,
                                                  minLength=self.minLength)
            fileForClustering = screenedFasta
        # Otherwise if masking is disabled, we'll use unique-ify and 
        #    pre-cluster our sequences
        else:
            uniqueFile, nameFile = self.uniqueSequences( filteredFile )
            preclusteredFile, nameFile = self.preclusterSequences( uniqueFile, nameFile )
            fileForClustering = preclusteredFile
        # If enabled, calculate sequence distances and cluster
        if self.enableClustering:
            distanceMatrix = self.calculateDistanceMatrix( fileForClustering )
            listFile = self.clusterSequences( distanceMatrix )
        # If enabled, generate a consensus for each cluster from above
        if self.enableConsensus:
            clusterListFile = self.separateClusterSequences( listFile, fastqFile )
            consensusFile = self.generateConsensusSequences( clusterListFile )
            self.cleanupConsensusFolder( consensusFile )
            selectedFile = self.selectFinalSequences( consensusFile )
"""

if __name__ == '__main__':
    rdnap = rDnaPipeline()
    rdnap()
