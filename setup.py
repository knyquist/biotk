#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, Extension, find_packages
import sys

if ("install" in sys.argv) and sys.version_info < (2, 7, 0):
    print "biochemistry-toolkit requires Python 2.7"
    sys.exit(-1)

globals = {}
execfile("biotk/__init__.py", globals)
__VERSION__ = globals["__VERSION__"]

setup(
    name='biotk',
    version=__VERSION__,
    author='Pacific Biosciences',
    author_email='devnet@pacificbiosciences.com',
    description="A Python library for doing custom biochemistry analysis",
    license='LICENSES.txt',
    packages=find_packages('.'),
    package_dir={'':'.'},
    zip_safe=False,
    scripts=['biotk/libs/BurstMetrics.py'],
    entry_points={'console_scripts':[
        'taulysis = biotk.tools.taulysis.main:main',
        'bammend = biotk.tools.bammend.main:main',
        'sequistory = biotk.scripts.screening_history:main',
        'pbisave = biotk.scripts.extractIndexMetrics:main',
        'curtail = biotk.scripts.rejectBasesByTime:main',
        'windchip = biotk.scripts.filterChips:main',
        'heatmap = biotk.scripts.plotHeatmap:main']})
