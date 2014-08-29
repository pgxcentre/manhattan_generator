#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip

import os
from setuptools import setup

# Check the Python version
import sys
major, minor, micro, s, tmp = sys.version_info
if major==2 and minor<7 or major<2:
    raise SystemExit("""Requires Python 2.7 or later.""")
if major==3:
    raise SystemExit("""Doesn't work on Python 3...""")

setup(name="manhattan_generator",
      version="1.6",
      description="Creation of beautiful Manhattan plots",
      author="Louis-Philippe Lemieux Perreault",
      author_email="louis-philippe.lemieux.perreault@statgen.org",
      url="http://www.statgen.org",
      license="GPL",
      scripts=[os.path.join("scripts", "manhattan_generator"),],
      install_requires=["matplotlib >=1.3.1", "numpy >= 1.8.0"],
      classifiers=['Operating System :: Linux',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 2.7'])
