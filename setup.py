#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.md') as history_file:
    history = history_file.read()

requirements = ['biopython', 'scipy','numpy','pandas','tables']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="VCFMIX allows identification of, and quantification of mixtures of, high quality bases within VCF files. It is designed for use with reference mapped bacterial sequences.",
    install_requires=requirements,
    license="GNU General Public License v3.0",
    long_description_content_type='text/markdown',
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='VCFMIX',
    name='VCFMIX',
    packages=find_packages(include=['VCFMIX', 'VCFMIX.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/AlexOrlek/VCXMIX',
    version='0.1.0',
    zip_safe=False,
)
