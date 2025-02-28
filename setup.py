# -*- coding: utf-8 -*-

"""setuptools control."""

import re
from setuptools import setup, find_packages
 
 
version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('tairdbsuite/tairdb.py').read(),
    re.M
    ).group(1)
 
 
with open("README.rst", "rb") as f:
    long_descr = f.read().decode("utf-8")
 
 
setup(
    name = "tairdbsuite",
    version = version,
    packages = find_packages(),
    description = "Python command line application for creating and querying a tair database.",
    long_description = long_descr,
    author = "Christian Goeschl",
    author_email = "christian.goeschl@gmi.oeaw.ac.at",
#     url = "http://gehrcke.de/2014/02/distributing-a-python-command-line-application",
    
    install_requires = [
        'numpy>=1.8.2',
        'scipy>=0.16.0',
        'pandas>=0.16.2',
        'sqlalchemy>=1.0.11',
        'statsmodels>=0.6.1'
    ],
      
    entry_points = {
        "console_scripts": ['tairdbsuite = tairdbsuite.tairdbsuite:main']
        },
    )
