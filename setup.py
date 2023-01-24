#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2022--, Panpiper development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="panpiper", 
    version="0.1.1",
    license='BSD-3-Clause',
    author="Renee Oles",
    author_email="roles@health.ucsd.edu",
    description="Panpiper: snakemake workflow for bacterial isolate analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rolesucsd/Panpiper",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Development Status :: 1 - Beta"],
    python_requires='>=3.5',
    install_requires=[
        "setuptools"],
    entry_points = {
        'console_scripts': ['panpiper = panpiper.main:cli']
    },
    include_package_data=True,
    zip_safe=False,
    package_data = {'panpiper': ['workflow/*']}
)
