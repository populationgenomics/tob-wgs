#!/usr/bin/env python

"""
Setup script for the pytob package
"""

import setuptools

with open('./README.md', 'r') as readme:
    long_description = readme.read()

setuptools.setup(
    name='pytob',
    packages='pytob',
    version='0.0.1',  # automatically updated by bump2version
    description='Python analysis modules for TOB-WGS.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url=f'https://github.com/populationgenomics/tob-wgs',
    license='MIT',
    author='Centre for Population Genomics',
    author_email='genomic-analysis-team@populationgenomics.org.au',
    include_package_data=True,
    zip_safe=False,
    keywords='bioinformatics',
    classifiers=[
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)
