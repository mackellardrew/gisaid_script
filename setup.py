#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = ['pytest>=3', ]

setup(
    author="Drew MacKellar",
    author_email='drew.mackellar@doh.wa.gov',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A package for formatting metadata for SARS-CoV-2 metadata",
    entry_points={
        'console_scripts': [
            'gisaid_script=gisaid_script.cli:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='gisaid_script',
    name='gisaid_script',
    packages=find_packages(include=['gisaid_script', 'gisaid_script.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/mackellardrew/gisaid_script',
    version='0.1.0',
    zip_safe=False,
)
