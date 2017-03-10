#!/usr/bin/env python

from setuptools import setup

import sys

setup(name="MICe-lab",
      version='0.1',
      #install_requires=['numpy', 'scipy', 'pyminc'],
      install_requires=[
        'ConfigArgParse>=0.11',
        'numpy',
        'pyminc',
        'pydpiper',
        'typing',
      ],
      scripts=['python/saddle_recon/saddle_recon.py'],
      )
