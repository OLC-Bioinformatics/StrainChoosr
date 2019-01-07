#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='diversitree',
    version='0.0.1',
    packages=find_packages(),
    scripts=['diversitree/diversitree.py'],
    author='Andrew Low',
    author_email='andrew.low@canada.ca',
    url='https://github.com/lowandrew/DiversiTree',
    install_requires=['biopython',
                      'scipy',
                      'matplotlib']
)
