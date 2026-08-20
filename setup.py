#!/usr/bin/env python
# -*- encoding: utf-8 -*-
"""Compatibility shim.

All packaging metadata now lives in pyproject.toml. This file only exists so
that legacy `python setup.py ...` invocations still find a setup() call.
"""
from setuptools import setup

setup()
