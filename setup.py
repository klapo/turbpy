import os
import re
import sys
import warnings
try:
    from setuptools import setup
    from setuptools.extension import Extension
except:
    from distutils.core import setup
    from distutils.extension import Extension

setup(name='turbpy',
      author="Karl Lapo",
      author_email="karl.lapo@uni-bayreuth.de",
      description="Turbulence related functions for land surface models",
      version='1.1',
      packages=['turbpy'],
      package_data={'turbpy': ['default_params.yml']},
      )
