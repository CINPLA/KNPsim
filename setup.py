# -*- coding: utf-8 -*-
from setuptools import setup
import os

from setuptools import find_packages


# long_description = open("README.md").read()

install_requires = []

setup(
    name="knpsim",
    packages=find_packages(),
    include_package_data=True,
    version="0.0.1"
)
