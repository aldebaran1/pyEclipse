#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on April 26 2021

@author: Sebastijan Mrak <smrak@bu.edu>
"""

from setuptools import setup


setup(name='pyEclipse',
      description='Software framework for computation of solar eclipses on the Earth.',
      author='Sebastijan Mrak',
      url='https://github.com/aldebaran1/pyEclipse.git',
      install_requires=['numpy', 'scipy', 'pyEphem', 'futures3', 'sunpy', 'zeep', 'lxml', 'drms', 'wget'],
      packages=['eclipse']
)