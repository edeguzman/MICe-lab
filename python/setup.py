#!/usr/bin/env python

from setuptools import setup, find_packages

import sys

setup(name="MICe-lab",
      version='0.18',
      description='Scripts and pipelines that are specific to the Mouse Imaging Centre lab.',
      long_description='Scripts and pipelines that are specific to the Mouse Imaging Centre lab.',
      url='https://github.com/Mouse-Imaging-Centre/MICe-lab',
      install_requires=[
        'ConfigArgParse>=0.11',
        'numpy',
        'pyminc',
        'typing',
        'scipy',
        'snakemake',
        'matplotlib',
        'Pillow'
      ],
      dependency_links=['https://github.com/opencv/opencv/archive/3.1.0.tar.gz',
	'https://github.com/Mouse-Imaging-Centre/pydpiper'],
      packages=find_packages(),
      scripts=['saddre/saddle_recon.py',
               'convert_CT_image.py',
               'Snakefile',
               'crop_to_brain.py',
               'mri_python/mri_recon.py',
               'mri_python/fse3dmice_recon.py',
               'tissvis/tissue_vision.py',
               'tissvis/TV_stitch.py',
               'tissvis/stacks_to_volume.py'
               ],
      )
