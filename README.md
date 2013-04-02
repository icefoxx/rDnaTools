# rDnaTools #

rDnaTools is a python package of tools and pipelines for working with
ribosomal DNA sequence data generated with the PacBio(R) SMRT sequencing.
rDnaTools works by wrapping existing tools from microbial ecology,
primarily the Mothur suite of utilities.

Currently rDnaTools implements a single pipeline for the export, filtering,
and cluster of 16S sequences.  Future releases will include automated
pipelines for other use-cases, as well as the capability for users to
script their own pipelines for rDNA sequence analysis.

Though primarily intended for use in analyzing 16S rDNA sequences, the
same tools and approaches should apply equally well to 18S, 23S, or ITS
sequences, provided that suitable reference sequences are supplied.

## Requirements ##

The core functionality of rDnaTools is built upon Python2.7 using the
pbcore framework for accessing PacBio data files.  In addition rDnaTools
wraps the functionality from a number of stand-alone commandline tools
that must available for the package to function
* Python 2.7 (www.python.org)
* pbcore (www.github.com/PacificBiosciences/pbcore)
* Blasr (www.github.com/PacificBiosciences/blasr)
* Mothur (www.mothur.org)

rDnaTools also includes the capability to generate high-quality consensus 
sequences from clusters of ribosomal DNA reads using algorithms 
included in the Pacific Biosciences SMRT Analysis suite.
* SMRTanalysis (www.pacbiodevnet.com)

## Primary Tools ##

rDnaTools contains two primary tools for analyzing rDNA sequence data
* rDnaPipeline
* rDnaResequencer

The rDnaPipeline takes as an input PacBio sequence data from ribosomal
DNA amplicons, in either FOFN, BAS.H5, FASTA, or FASTQ formate, and runs a
sequential series of analyses, similar to Mothur`s Batch Mode.

The rDnaResequencer tool generates consensus sequences from clusters of
PacBio ribosomal DNA sequences.  The Resequencer tool can be called either
independently, or as part of the rDnaPipeline.

## Citation ##

rDnaTools would not have been possible were it not for the hard work of the
existing Microbial Ecology community, and their existing tools for analyzing
ribosomal DNA sequence data.  Since the core of the analyses wrapped by 
rDnaTools come from the Mothur suite, please cite their publication if you 
use rDnaTools in your work:

Schloss, P.D., et al., Introducing mothur: Open-source, platform-independent, 
community-supported software for describing and comparing microbial communities. 
Appl Environ Microbiol, 2009. 75(23):7537-41.
