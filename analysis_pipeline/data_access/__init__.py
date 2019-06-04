#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides access to supernova light-curve data from DES, SDSS,
and CSP. Data is downloaded automatically if it is not locally available,
including filter transmission curves. This package will temporarily register
filter transmission curves with SNCosmo using the naming scheme
`91bg_proj_<survey name>_<filter name>`."""

from . import csp, des, sdss