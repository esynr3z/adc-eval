#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Analog <-> Digital converters behavioral models
"""

import numpy as np


def analog2digital(sig_f, sample_freq=1e6, sample_n=1024, sample_bits=8, vref=3.3, noisy_lsb=1):
    sample_quants = 2 ** sample_bits
    sample_prd = 1 / sample_freq
    t = np.arange(0, sample_n * sample_prd, sample_prd)
    dv = vref / sample_quants
    samples = np.rint(sig_f(t) / dv).astype(int)
    # apply noise
    if noisy_lsb:
        noise = np.random.normal(0, 2 ** (noisy_lsb - 1), size=sample_n).astype(int)
        samples += noise
        # to be sure that samples fit the range:
        samples[np.argwhere(samples >= sample_quants)] = sample_quants - 1
        samples[np.argwhere(samples < 0)] = 0
    return samples


def digital2analog(samples, sample_bits=8, vref=3.3):
    quants = 2 ** sample_bits
    dv = vref / quants
    return samples * dv
