#!/usr/bin/env python

import logging
import subprocess

from .._utils import which, fileExists, isExe

__version__ = "0.1"

VALID_COMMANDS = frozenset(['fastq.info', 'align.seqs', 'chimera.uchime',
                            'screen.seqs', 'remove.seqs', 'filter.seqs', 
                            'summary.seqs', 'dist.seqs', 'cluster', 
                            'set.logfile'])

NUMPROC_COMMANDS = frozenset(['align.seqs', 'chimera.uchime', 'screen.seqs',
                              'filter.seqs', 'dist.seqs'])

VALID_PARAMS = frozenset(['fasta', 'fastq', 'qfile', 'reference', 'name',
                          'flip', 'start', 'end', 'minlength', 'processors',
                          'vertical', 'calc', 'output', 'phylip', 'method',
                          'accnos'])

class MothurCommand(object):
    """
    A class for representing individual Mothur Commands
    """
    def __init__(self, command, parameters):
        # Validate the selected Mothur Command
        try:
            assert command.lower() in VALID_COMMANDS
            self._command = command.lower()
        except AssertionError:
            ValueError("Invalid Mothur command")
        # Validate the arguments for the command
        self._parameters = {}
        for param, value in parameters.iteritems():
            try: 
                assert param in VALID_PARAMS
                self.parameters[param] = value
            except AssertionError:
                ValueError('Invalid Mothur argument')

    @property
    def command(self):
        return self._command

    @property
    def parameters(self):
        return self._parameters

    @property
    def parameterString(self):
        return ", ".join(["%s=%s"%(p,v) for p, v in self.parameters.iteritems()
                                        if v is not None])

    def __str__(self):
        return "%s(%s)" % (self.command, 
                           self.parameterString)



class MothurJob(object):
    """
    A class for representing Mothur calls of 1-or-more commands 
    """
    def __init__(self, mothurExe, commands, stdout=None, stderr=None):
        # Validate the supplied Mothur executable
        try:
            assert isExe(mothurExe)
            self._mothur = mothurExe
        except AssertionError:
            raise ValueError("Mothur executable is not valid!")
        # Validate the supplied Mothur commands
        self._commands = []
        for command in commands:
            try:
                assert isinstance(command, MothurCommand)
                self._commands.append( command )
            except:
                raise ValueError("Argument is not a MothurCommand!")
        self.stdout = stdout
        self.stderr = stderr

    @property
    def mothur(self):
        return self._mothur

    @property
    def commands(self):
        return self._commands

    @property
    def commandString(self):
        return "; ".join( map(str, self.commands) )

    @property
    def callString(self):
        return '#%s' % self.commandString

    def __str__(self):
        return '%s "%s"' % (self.mothur,
                            self.callString)

    def __call__(self):
        logging.info('Executing the following command:\n\t%s' % str(self))
        p = subprocess.Popen( [self.mothur, self.callString], 
                              stdout=self.stdout, stderr=self.stderr )
        stdout, stderr = p.communicate()
        logging.info('Mothur command exited successfully')



class MothurFactory(object):
    """
    A Factory-style tool for creating run-able MothurCommand objects
    """
    #####################
    # Variable Defaults #
    #####################

    NUM_PROC = 1

    ##########################
    # Initialization Methods #
    ##########################

    def __init__(self, mothurExe=None, numProc=None, stdout=None, stderr=None):
        if mothurExe is None:
            self.initializeFromArgs()
        else:
            self.initializeFromCall(mothurExe, numProc, stdout, stderr)
        self.validateSettings()

    def initializeFromArgs(self):
        import argparse
        desc = 'A tool for calling Mothur commands from Python'
        parser = argparse.ArgumentParser(description=desc)
        parser.add_argument('-n', '--num_processes', type=int,
                            metavar='INT', dest='numProc',
                            default=self.NUM_PROC,
                            help="Number of processors to use")
        parser.add_argument('--mothur', metavar='MOTHUR_PATH', 
                            dest='mothurExe',
                            help="Specify the path to the Mothur executable")
        parser.add_argument('--stdout', metavar='OUT',
                            type=argparse.Filetype('a'),
                            help="Pipe Mothur's STDOUT to specified file")
        parser.add_argument('--stderr', metavar='ERR',
                            type=argparse.FileType('a'),
                            help="Pipe Mothur's STDERR to specified file")
        parser.add_argument('--debug', action='store_true',
                            help="Turn on DEBUG message logging")
        parser.add_argument('--debug', action='store_true',
                            help="Turn on DEBUG message logging")
        args = parser.parse_args()
        self.__dict__.update( vars(args) )

    def initializeFromCall(self, mothurExe, numProc, stdout, stderr):
        self.mothurExe = mothurExe
        self.numProc = numProc
        # Set the output handle for STDOUT
        if stdout is None:
            self.stdout = None
        else:
            self.stdout = open(stdout, 'a')
        # Set the output handle for STDERR
        if stderr is None:
            self.stderr = None
        else:
            self.stderr = open(stderr, 'a')

    def validateSettings(self):
        # Search for Mothur executable if not supplied
        if self.mothurExe is None:
            self.mothur = which('mothur')
            if self.mothur is None:
                raise OSError('Mothur executable not found!')
        # If an argument was given, check that it's executable
        elif isExe( self.mothurExe ):
            self.mothur = self.mothurExe
        else:
            raise OSError('Supplied Mothur executable not valid!')
        # Validate the Num_Processes argument
        try:
            assert self.numProc >= 1
        except AssertionError:
            raise ValueError("Number of processes must be >= 1!")
        self.numProc = str(self.numProc)

    ####################
    # Instance Methods #
    ####################

    def createJob(self, command, parameters, logFile=None):
        commands = []
        logging.info('Creating a MothurCommand for "%s"' % command)
        # Check the logFile and create it if needed
        if logFile is None:
            logging.info('No log-file specified for this job')
        else:
            logging.info('Setting the log-file to "%s"' % logFile)
            logParams = {'name':logFile}
            logCommand = MothurCommand('set.logfile', logParams)
            commands.append( logCommand )
        parameters = self.addDefaultParameters( command, parameters )
        mainCommand = MothurCommand( command, parameters )
        commands.append( mainCommand )
        return MothurJob(self.mothur, commands, self.stdout, self.stderr)

    def addDefaultParameters(self, command, parameters):
        logging.info('Checking default parameters for "%s"' % command)
        if command in NUMPROC_COMMANDS and 'processors' not in parameters:
            logging.info('Adding default value from "numProc" for "processors"')
            parameters['processors'] = str(self.numProc)
        return parameters

    def runJob(self, command, parameters, logFile=None):
        job = self.createJob(command, parameters, logFile)
        job()

if __name__ == '__main__':
    mcm = MothurFactory()
    print "MothurFactory Created"
