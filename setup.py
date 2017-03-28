#!/usr/bin/env python

from setuptools import setup

import sys

setup(name="MICe-lab",
      version='0.1',
      #install_requires=['numpy', 'scipy', 'pyminc'],
      scripts=[
               "python/convert_CT_image.py",
               ],
      )
