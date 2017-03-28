#!/usr/bin/env python

from setuptools import setup

import sys

setup(name="MICe-lab",
      version='0.1',
      install_requires=['snakemake'], #'numpy', 'scipy', 'pyminc'],
      scripts=[
               "python/convert_CT_image.py",
               "python/Snakefile"
               ],
      )
