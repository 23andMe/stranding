#!/usr/bin/env python
import sys
from setuptools import setup, find_packages

VERSION = '0.1.8'

PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3

if PY2:
    readme = open('README.md').read()
    license = open('LICENSE').read()
elif PY3:
    readme = open('README.md', encoding='utf-8').read()
    license = open('LICENSE', encoding='utf-8').read()

setup(
    name='stranding',
    version=VERSION,
    url='https://github.com/23andMe/stranding',
    download_url = 'https://github.com/23andMe/stranding/tarball/%s' % VERSION,
    author='23andMe Engineering',
    author_email=['mstrand@23andme.com'],
    description='',
    long_description=readme,
    entry_points={'console_scripts': []},
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['seqseek', 'biopython'],
    classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ]
)
