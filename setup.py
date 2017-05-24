#!/usr/bin/env python
import os
import sys
import re

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

# Handle encoding
major, minor1, minor2, release, serial = sys.version_info
if major >= 3:
    def rd(filename):
        f = open(filename, encoding="utf-8")
        r = f.read()
        f.close()
        return r
else:
    def rd(filename):
        f = open(filename)
        r = f.read()
        f.close()
        return r

setup(
    name='moscatel',
    packages =['moscatel'],
    version="0.1.1",
    author='John Livingston and Jerome de Leon',
    author_email = 'jpdeleon.bsap@gmail.com',
    url = 'https://github.com/jpdeleon/moscatel',
    license = ['GNU GPLv3'],
    description ='A simple pipeline for obtaining light curve of MUSCAT data.',
    long_description=rd("README.md") + "\n\n"
                    + "---------\n\n",
    package_dir={"moscatel": "moscatel"},
    scripts=['scripts/moscatel-phot', 'scripts/moscatel-analysis', 'scripts/moscatel-init'],
    include_package_data=True,
    keywords=['muscat','multi-color','transit photometry'],
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python'
        ],
    install_requires = ['numpy', 'pandas', 'matplotlib', 'astropy', 'photutils', 'tqdm', 'scikit-image'],
)
