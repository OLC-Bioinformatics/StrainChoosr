#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='strainchoosr',
    version='0.0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'strainchoosr = strainchoosr.strainchoosr:main',
        ],
    },
    author='Andrew Low',
    author_email='andrew.low@canada.ca',
    url='https://github.com/lowandrew/StrainChoosr',
    install_requires=['scipy',
                      'pytest',
                      'ete3']
)
