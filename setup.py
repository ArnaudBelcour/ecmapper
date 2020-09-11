#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['biopython',
                'cobra']

setup(
    author="Arnaud Belcour",
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="test",
    entry_points={
        'console_scripts': [
            'ecmapper=ecmapper.__main__:cli',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    include_package_data=True,
    keywords='ecmapper',
    name='ecmapper',
    packages=find_packages(include=['ecmapper', 'ecmapper.*']),
    zip_safe=False,
)