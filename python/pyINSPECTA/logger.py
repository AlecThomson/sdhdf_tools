#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Logging module for pyINSPECTA"""

import logging

# Create logger
logger = logging.getLogger('pyINSPECTA')
logger.setLevel(logging.DEBUG)

# Create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# Create formatter
formatter = logging.Formatter('%(name)s - %(levelname)s: %(message)s')

# Add formatter to ch
ch.setFormatter(formatter)

# Add ch to logger
logger.addHandler(ch)
