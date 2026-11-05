#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name='strainchoosr',
    version='0.1.7',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'strainchoosr = strainchoosr.strainchoosr:main',
            'strainchoosr_gui = strainchoosr.strainchoosr_gui:main',
            'strainchoosr_drawimage = strainchoosr.strainchoosr_gui:draw_image_wrapper'
        ],
    },
    author='Adam Koziol',
    author_email='adam.koziol@inspection.gc.ca',
    url='https://github.com/OLC-Bioinformatics/StrainChoosr',
    python_requires='>=3.5',
    install_requires=[
        'ete3>=3.1.3',
        'numpy>=1.15',
        'six>=1.16.0',
        'legacy-cgi>=2.6.1; python_version >= "3.13"',
        'importlib-metadata>=1.0; python_version < "3.8"',
    ],
    extras_require={
        'test': ['pytest'],
        'visualization': [
            'PyQt5>=5.15.11; python_version >= "3.8"',
            'PyQt5>=5.11.3,<5.12; python_version < "3.8"',
        ],
    },
)
