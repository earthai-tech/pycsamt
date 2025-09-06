#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 10 18:49:02 2022

@author: daniel
"""
import os 
import tempfile
# set path for demo and testing module
# Get the repository dir like https://github.com/WEgeophysics/pycsamt /
#C:github.com/WEgeophysics/pyCSAMT  
PYCSAMT_ROOT = os.path.normpath(
    os.path.abspath(
        os.path.dirname(
            os.path.dirname(__file__)
        )
    )
)

EDI_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/edi'))
AVG_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/avg'))
J_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/j'))
DRILL_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/drill_examples_files'))
OCCAM2D_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/occam2d'))
STN_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/stn_profiles'))
CONFIG_DATA_DIR = os.path.normpath(
    os.path.join(PYCSAMT_ROOT , 'data/_conffiles'))

SYSTEM_TEMP_DIR = tempfile.gettempdir()
NEW_TEMP_DIR=tempfile.mkdtemp(prefix="pycsamt_tmpdir_")

