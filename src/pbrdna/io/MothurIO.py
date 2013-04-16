import sys
import csv

from collections import namedtuple, Counter

SummaryRecord = namedtuple('SummaryRecord', 'seqname start end nbases ambigs polymer numSeqs')

class SummaryReader(object):
    """
    A tool for parsing data from Mothur Summary files
    """
    DEFAULT_FRAC = 0.8

    def __init__(self, summary=None, fraction=None):
        if summary is None:
            self.initializeFromArgs()
        else:
            self.initializeFromCall(summary, fraction)
        self.validateSettings()
        self.parseSummaryData()

    def initializeFromArgs(self):
        import argparse
        desc = "A tool for parsing data from Mothur Summay Files"
        parser = argparse.ArgumentParser(description=desc)
        parser.add_argument('summary', metavar='FILE',
                            help='Summary file to parse')
        parser.add_argument('-f', '--fraction', metavar='FLOAT',
                            type=float, default=self.DEFAULT_FRAC,
                            help='Minimum fraction of the full read to ' + \
                                 'calculate f')
        args = parser.parse_args()
        self.__dict__.update( vars(args) )

    def initializeFromCall(self, summary, fraction):
        self.summary = summary
        if fraction is None:
            self.fraction = self.DEFAULT_FRAC
        else:
            self.fraction = fraction

    def validateSettings(self):
        # Validate the Summary file
        try:
            assert self.summary.endswith('.summary')
        except AssertionError:
            raise ValueError('"%s" is not a valid summary file!' % self.summary)
        # Validate the fraction argument
        try:
            assert 0.1 < self.fraction <= 1.0
        except AssertionError:
            raise ValueError('Fraction value must be >0.1 and <=1.0!')
        self.start = None
        self.end = None

    def parseSummaryData(self):
        self.summaryData = []
        with open(self.summary, 'r') as handle:
            for record in map(SummaryRecord._make, csv.reader(handle, delimiter='\t')):
                self.summaryData.append(record)

    def parseStart(self):
        counter = Counter( [record.start for record in self.summaryData] )
        return int(counter.most_common(1)[0][0])

    def parseEnd(self):
        counter = Counter( [record.end for record in self.summaryData] )
        return int(counter.most_common(1)[0][0])

    def getFullLengthPositions(self):
        if self.start is None:
            self.start = self.parseStart()
        if self.end is None:
            self.end = self.parseEnd()
        return (self.start, self.end)

    def getAllowedPositions(self):
        if self.start is None or self.end is None:
            self.getFullLengthPositions()
        length = self.end - self.start
        margin = (1-self.fraction) / 2
        minimumEnd = int(self.end - (margin * length))
        maximumStart = int(self.start + (margin * length))
        return (maximumStart, minimumEnd)

if __name__ == '__main__':
    parser = SummaryReader()
    print parser.getFullLengthPositions()
