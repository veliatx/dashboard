# -*- coding: utf-8 -*-

from os.path import abspath, dirname
from sys import path
from setuptools import setup, find_packages

setup(
    name='dashboard',
    version='0.0.1',
    description="""Velia Dashboard facilitates visualization and exploration of sORF related data.""",
    author='Rob Foreman',
    author_email='rob@veliatx.com',
    classifiers=[
        'Programming Language :: Python :: 3.9',
    ],
    keywords='microproteins',
    packages=find_packages(),
    install_requires=[],
    entry_points = {
        'console_scripts': ['dashboard_update=dashboard.etl.update_cache:update_cache', 
                            'run_protein_search_tools=dashboard.etl.run_protein_search_tools:run_protein_search_tools'],
    }
)