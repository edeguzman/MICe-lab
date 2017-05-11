#!/usr/bin/env python

from setuptools import setup

import sys

setup(name="MICe-lab",
      version='0.12',
      install_requires=[
        'ConfigArgParse>=0.11',
        'numpy',
        'pyminc',
        'typing',
        'snakemake'
      ],
      packages=['python', 'python.saddle_recon'],
      scripts=['python/saddle_recon/saddle_recon.py',
               'python/convert_CT_image.py',
               'python/Snakefile'
               ],
      )
