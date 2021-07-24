#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Some basic signal functions
"""

import numpy as np


def sin(t, peak=1.5, offset=1.65, freq=1e3, ph0=0):
    return offset + peak * np.sin(ph0 + 2 * np.pi * freq * t)


def noise(t, mean=0, std=0.1):
    return np.random.normal(mean, std, size=len(t))
