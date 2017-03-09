#!/usr/bin/env python
import io
from setuptools import setup, find_packages

VERSION = '0.1.9'

with io.open('README.md', encoding='utf-8') as f:
    readme = f.read()
with io.open('LICENSE', encoding='utf-8') as f:
    license = f.read()


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
