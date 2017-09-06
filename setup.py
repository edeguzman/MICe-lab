#!/usr/bin/env python

from setuptools import setup

import sys

setup(name="MICe-lab",
      version='0.15',
      install_requires=[
        'ConfigArgParse>=0.11',
        'numpy',
        'pyminc',
        'typing',
        'snakemake',
        'matplotlib'
      ],
      packages=['python', 'python.saddle_recon','python.mri_python'],
      scripts=['python/saddle_recon/saddle_recon.py',
               'python/convert_CT_image.py',
               'python/Snakefile',
               'python/crop_to_brain.py',
               'python/mri_python/mri_recon.py',
               ],
      )
