#! /usr/bin/env python

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

from pbtools.pbrdna.io.FastqIO import FastqReader, FastqWriter
from pbtools.pbrdna.io.FastaIO import FastaRecord, FastaWriter  
from pbtools.pbrdna._utils import createDirectory, validateInputFile, validateFloat
from pbtools.pbrdna._utils import which, getZmw

DEFAULT_DIST = 0.03
MIN_FULL_LENGTH = 1400

class ClusterSeparator(object):
    """
    A tool for resequencing clusters of rDNA sequences    
    """

    ##########################
    # Initialization Methods #
    ##########################

    def __init__(self, listFile=None, ccsFile=None, distance=None):
        if listFile is None or ccsFile is None:
            self.initializeFromArgs()
        else:
            self.initializeFromCall(listFile, ccsFile, distance)
        self.validateSettings()

    def initializeFromArgs(self):
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('listFile', metavar='FILE',
                            help="Mothur list file of cluster data")
        parser.add_argument('ccsFile', metavar='FILE',
                            help="Fasta or Fastq file of CCS sequences")
        parser.add_argument('-d', '--distance', type=float,
                            default=DEFAULT_DIST, metavar='FLOAT',
                            help="Distance at which to cluster sequences")
        parser.add_argument('-m', '--minRefLength', type=int,
                            default=MIN_FULL_LENGTH, metavar='INT',
                            help="Minimum length to allow for reference sequences")
        parser.add_argument('--outputFastq', action='store_true',
                            help="Output FASTQ files for each cluster")
        parser.add_argument('--outputReference', action='store_true',
                            help="Output reference sequences for each cluster")
        parser.add_argument('-o', '--outputDir', default='reseq',
                            help="Specify a directory for output files")
        args = parser.parse_args()
        self.__dict__.update( vars(args) )

    def initializeFromCall(self, listFile, seqFile, distance):
        self.listFile = listFile
        self.sequenceFile = seqFile
        self.distance = DEFAULT_DIST if distance is None else distance
        self.outputDir = 'reseq'

    def validateSettings(self):
        # Check the values of the supplied input files
        self.listFile = validateInputFile( self.listFile, ['.list'])
        self.ccsFile = validateInputFile( self.ccsFile, ['.fq', '.fastq'])
        # Check the value of the supplied distance
        self.distance = validateFloat( self.distance, 
                                       minValue=0.001, 
                                       maxValue=0.5)
        # Create the output directory if needed and move into it
        createDirectory( self.outputDir )
        os.chdir( self.outputDir )

    #################
    # Class Methods #
    #################

    @classmethod
    def convertDistance(cls, distance):
        try:
            distance = 'unique' if distance == 'unique' else float(distance)
        except:
            raise ValueError('"%s" is not a valid distance!' % parts[0])
        return distance

    ####################
    # Instance Methods #
    ####################

    def parseSequenceData(self):
        self.sequenceData = {}
        for fastqRecord in FastqReader( self.ccsFile ):
            zmw = getZmw( fastqRecord.name )
            self.sequenceData[zmw] = fastqRecord

    def parseDistances(self):
        distances = []
        with open( self.listFile, 'r' ) as handle:
            for line in handle:
                parts = line.split()
                distance = self.convertDistance( parts[0] )
                distances.append( distance )
        return distances

    def selectDistance(self, distances):
        # If our selected distance is present, simply return it
        if self.distance in distances:
            return self.distance
        # Otherwise find the largest clustering distance smaller than 
        #    the specified distance and return that
        possible = [d for d in distances if d != 'unique']
        smaller = [d for d in possible if d < self.distance]
        if not smaller:
            raise ValueError('No valid clustering distances found!')
        return max(smaller)

    def parseClusters( self, targetDist ):
        with open( self.listFile, 'r' ) as handle:
            for line in handle:
                # Skip lines until we find the target distance
                parts = line.split()
                currDist = self.convertDistance( parts[0] )
                if currDist != targetDist:
                    continue
                # Check that the number of clusters is concordant
                clusterCount = int(parts[1])
                clusters = parts[2:]
                assert len(clusters) == clusterCount
                # Convert the strings of clusters to Lists and return
                clusters = [c.split(',') for c in clusters]
                return clusters

    def trimClusterNames(self, clusters):
        trimmed = []
        for cluster in clusters:
            cluster = [getZmw(c) for c in cluster]
            trimmed.append( frozenset(cluster) )
        return trimmed

    def getClusterReads(self, cluster):
        reads = []
        for ccsZmw in cluster:
            try:
                ccsRead = self.sequenceData[ccsZmw]
            except KeyError:
                #raise Warning("No CCS read found for '%s', skipping..." % ccsZmw)
                continue
            reads.append( ccsRead )
        return reads

    def outputClusterFastq(self, reads, count):
        fastqFile = 'cluster%s.fastq' % count
        if os.path.exists( fastqFile ):
            return fastqFile
        # Rename the "Reference" sequence to the cluster
        with FastqWriter( fastqFile ) as handle:
            for fastqRecord in reads:
                handle.writeRecord( fastqRecord )
        return fastqFile

    def outputClusterFasta(self, reads, count):
        fastaFile = 'cluster%s.fasta' % count
        if os.path.exists( fastaFile ):
            return fastaFile
        # Rename the "Reference" sequence to the cluster
        with FastaWriter( fastaFile ) as handle:
            for fastqRecord in reads:
                fastaRecord = FastaRecord( fastqRecord.name,
                                           fastqRecord.sequence )
                handle.writeRecord( fastaRecord )
        return fastaFile

    def pickReference(self, reads):
        longReads = [read for read in reads
                          if len(read.sequence) > self.minRefLength]
        if len(longReads) == 1:
            return longReads[0]
        elif longReads:
            return self.findLowestErrorRead( reads )
        # If no 'full-length' reads are present, simply return the longest
        else:
            return self.findLongestRead( reads )

    def findLowestErrorRead( reads ):
        qvTuples = []

    def findLongestRead(self, reads):
        lengths = [len(read.sequence) for read in reads]
        maxLength = max(lengths)
        longestReads = [read for read in reads
                             if len(read.sequence) == maxLength]
        return longestReads[0]

    def outputClusterReference(self, reference, count):
        print "Creating reference sequence for Cluster #%s" % count
        referenceFile = 'cluster%s_reference.fasta' % count
        if os.path.exists( referenceFile ):
            return referenceFile
        # Rename the "Reference" sequence to the cluster
        referenceFasta = FastaRecord("Cluster%s" % count,
                                     reference.sequence)
        with FastaWriter( referenceFile ) as handle:
            handle.writeRecord( referenceFasta )
        return referenceFile

    def outputRepresentativeRead(self, representativeRead, count):
        print "Creating representative sequence file Cluster #%s" % count
        representativeFile = 'cluster%s_represent.fastq' % count
        if os.path.exists( representativeFile ):
            return representativeFile
        with FastqWriter( representativeFile ) as handle:
            handle.writeRecord( representativeRead )
        return representativeFile

    def __call__(self):
        self.parseSequenceData()
        # Select the appropriate distance, and parse the matching clusters
        distances = self.parseDistances()
        distance = self.selectDistance( distances )
        clusters = self.parseClusters( distance )
        # Trim the cluster neams and iterate, outputing each subset
        trimmedClusters = self.trimClusterNames( clusters )
        for count, cluster in enumerate( trimmedClusters ):
            count = str(count+1).zfill(4)
            print "Analyzing cluster #%s now..." % (count)
            reads = self.getClusterReads( cluster )
            if self.outputFastq:
                self.outputClusterFastq( reads, count )
            if self.outputReference:
                reference = self.pickReference( reads )
                self.outputReferenceFasta( reference, count )
            self.outputClusterFasta( reads, count )
        # Finally we combine and trim all of the output Files

if __name__ == '__main__':
    separator = ClusterSeparator()
    separator()
