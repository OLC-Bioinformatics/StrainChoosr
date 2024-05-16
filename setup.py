#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='strainchoosr',
    version='0.1.3',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'strainchoosr = strainchoosr.strainchoosr:main',
            'strainchoosr_gui = strainchoosr.strainchoosr_gui:main',
            'strainchoosr_drawimage = strainchoosr.strainchoosr_gui:draw_image_wrapper'
        ],
    },
    author='Andrew Low',
    author_email='andrew.low@canada.ca',
    url='https://github.com/lowandrew/StrainChoosr',
    #setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    install_requires=['pytest',
                      'ete3',
                      'PyQt5==5.11.3',  # Apparently required by ete3 for visualisation, but not listed in that setup.py...
                      ]
)
