#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SHDF Error classes"""

class VerificationError(Exception):
    """ Error raised if Verification fails """
    def __init__(self, msg):
        super().__init__(self, msg)