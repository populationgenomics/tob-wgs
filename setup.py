#!/usr/bin/env python

"""
Setup script
"""

import setuptools

PKG = 'tob-wgs'

setuptools.setup(
    name=PKG,
    packages=[PKG],
    version='0.0.1',  # automatically updated by bump2version
    description=f'Python analysis modules for {PKG}.',
    long_description=open('./README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/populationgenomics/{PKG}',
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
