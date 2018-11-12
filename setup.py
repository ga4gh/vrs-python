# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

short_description = "demonstrations for using the VMC data models"
long_description = open("README.rst").read()


setup(
    author = "Reece Hart",
    author_email = "reecehart@gmail.com",
    description = short_description.replace("\n", " "),
    license = "Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0)",
    long_description = long_description,
    name = "vmc",
    packages = find_packages(),
    package_data = {"vmc": ["_data/schema/*.json"]},
    url = "https://github.com/reece/vmc-python",
    use_scm_version = True,
    zip_safe = True,

    classifiers = [
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Topic :: Database :: Front-Ends",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
    ],

    keywords = [
        'bioinformatics',
    ],

    extras_require = {
        'dev': [
            "pytest",
            "pytest-cov",
            "restview",
            "tox",
            "yapf",
        ],
        'notebooks': [
            "ipython",
            "jupyter",
        ],
    },

    install_requires = [
        "biocommons.seqrepo>=0.4.2",
        "bioutils>=0.4rc",
        "hgvs>=1.2",
        "python-jsonschema-objects<0.3",
    ],

    setup_requires = [
        "pytest-runner",
        "setuptools_scm",
        "wheel",
    ],

    tests_require = [
        "pytest",
        "pytest-cov",
        "pyyaml",
    ],
)

# <LICENSE>
# Copyright 2016 Source Code Committers
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# </LICENSE>
