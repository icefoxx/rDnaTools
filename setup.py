from setuptools import setup, Extension, find_packages
import sys

desc = 'Tools for analyzing SMRT sequencing data from ribosomal ' + \
       'DNA amplicons (16S, 23S, ITS)'

if ("install" in sys.argv) and sys.version_info < (2, 7, 0):
    raise SystemExit("rDnaTools requires Python 2.7")

setup(
    name = 'rDnaTools',
    version='0.1.0',
    author='Brett Bowman',
    author_email='bbowman@pacificbiosciences.com',
    url='https://github.com/bnbowman/rDnaTools',
    description=desc,
    license=open('LICENSES.txt').read(),
    packages = find_packages('src'),
    package_dir = {'':'src'},
    scripts=['src/rDnaPipeline.py'],
    zip_safe = False,
    install_requires=[
        'h5py >= 2.0.1',
        'numpy >= 1.6.0',
        'pbcore >= 0.6.0',
        'pbtools.pbdagcon >= 0.2.1'
        ]
    )
