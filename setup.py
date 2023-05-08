#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()


with open('requirements.txt') as reqs_file:
    requirements = reqs_file.read().strip().split("\n")


test_requirements = ['pytest>=3', ]

setup(
    author="Matthew Turk",
    author_email='matthewturk@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A frontend for libyt",
    install_requires=requirements,
    license="BSD license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='yt_libyt',
    name='yt_libyt',
    packages=find_packages(include=['yt_libyt', 'yt_libyt.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/matthewturk/yt_libyt',
    version='0.1.0',
    zip_safe=False,
)
