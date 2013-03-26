Installing rDnaTools and Dependencies
=====================================

:Authors: Brett Bowman

:Version: 0.1

:Date: March 25th, 2013


Prerequisites:

* A UNIX-derived operating system
* gcc >4.4.3
* python 2.7
* git ( http://git-scm.com/ for installing some dependencies from github )
* blasr ( https://github.com/PacificBiosciences/blasr )

Optional:

* Quiver from SMRTAnalysis > v1.4 ( from http://pacbiodevnet.com ) 
  This is only required for creating cluster consensus sequences

Creating the Virtual Environment
--------------------------------

First we need to create a directory to store the virtual environment, and
export the path to a variable we can call later.  I prefer to keep my 
environments in a separate folder, i.e. ~/environments/RDNA_ENV, but
any clean directory with write-access should work::
    $ mkdir /some/path/to/your/RDNA_ENV
    $ export RDNA_HOME=/some/path/to/your/RDNA_ENV

Now we can create a clean Virtual Environment in the new directory, making
sure that we specify Python 2.7, and move into the directory::
    $ virtualenv -p /usr/bin/python2.7 $RDNA_HOME
    $ cd $RDNA_HOME

First we need to activate the new environment, so that when we install 
prerequisites they are applied only to the new Virtual Environment::
    $ source bin/activate

Installing Mothur
-----------------

Mothur is a widely used interactive tool for carrying out the analysis of
16S and other ribosomal DNA sequences.  For the parts of the 16S analysis
pipeline that don't require special tools to accomodate SMRT sequencing
data, rDnaTools simply wraps that functionality from Mothur.

To download and install local binaries in the Virtual Environment and
add them to the local PATH::
    $ cd bin/
    $ wget http://www.mothur.org/w/images/8/88/Mothur.cen_64.zip
    $ unzip Mothur.cen_64.zip -d Mothur_Dist
    $ cd Mothur_Dist/mothur
    $ cp -r * ../../
    $ cd $RDNA_HOME

Installing Python Library Dependencies
--------------------------------------

Next we need to install the public Python packages that the PacBio-specific 
packages require, of which there are two - ``numpy`` and ``h5py``.  If you 
know that current versions of these libraries are already installed in your
python environment, you should skip down to the next section.

If you don't have an installation of ``numpy`` ( http://www.numpy.org ) in
your current python environemtn, you can install it and it's dependencies 
with ``pip``::
    $ pip install numpy==1.6.2

Next we need ``h5py``, which in turn depends on ``libhdf5`` 
( http://www.hdfgroup.org/ftp/HDF5/prev-releases/ ).  In order to download 
and install ``libhdf5`` in the virtualenv::
    $ cd bin/
    $ wget http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.gz
    $ tar zxvf hdf5-1.8.9.tar.gz
    $ cd hdf5-1.8.9
    $ ./configure --prefix=$RDNA_HOME --enable-cxx
    $ make install
    $ cd $RDNA_HOME

Finally, we can install ``h5py`` ( http://h5py.googlecode.com/ ) to our
local python environment::
    $ cd bin/
    $ wget http://h5py.googlecode.com/files/h5py-2.0.1.tar.gz
    $ tar zxvf h5py-2.0.1.tar.gz
    $ cd zxvf h5py-2.0.1
    $ python setup.py build --hdf5=$RDNA_HOME
    $ python setup.py install
    $ cd $RDNA_HOME

Installing PacBio-specific Python Libraries
-------------------------------------------

Now we can finally install the PacBio python libraries 
and rDnaTools::
    $ pip install git+https://github.com/PacificBiosciences/pbcore.git
    $ pip install git+https://github.com/bnbowman/rDnaTools.git

If you do a ``pip freeze`` at this junction to inspect your installed libraries,
this is what you should see::
    $ pip freeze
    h5py==2.0.1
    numpy==1.6.2
    pbcore==0.6.0
    wsgiref=0.1.2

Install Other rDnaTools Prerequisites
-------------------------------------

The last requirement for rDnaTools is ``BLASR`` for the accurate mapping of
rDNA reads to their references.  The ``BLASR`` executable is included in the 
SMRT Analysis suite provided by PacBio, and  is also available on github. 
If you've already have a version of SMRT Analysis installed, you can make a 
local copy with the following command::
    $ cp $(which blasr) $RDNA_HOME/bin

For the final, optional step in the rDNA analysis pipeline is creating a high-
quality consensus sequence for each cluster using ``Quiver``, for which we need
a full installation of the SMRT Analysis installation.  Full installation
binaries and instructions can be found on PacBio's DevNet 
( http://pacbiodevnet.com )
