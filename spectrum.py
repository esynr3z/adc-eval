#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Some basic spectral analysis.

References:
  - Analog Devices MT-003 TUTORIAL
    "Understand SINAD, ENOB, SNR, THD, THD + N, and SFDR so You Don't Get Lost in the Noise Floor"
  - National Instruments Application Note 041
    "The Fundamentals of FFT-Based Signal Analysis and Measurement"
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def db2amp(db):
    """Decibels to amplitutde ratio"""
    return 10 ** (0.05 * db)


def amp2db(a):
    """Amplitutde ratio to decibels"""
    return 20 * np.log10(a)


def db2pow(db):
    """Decibels to power ratio"""
    return 10 ** (0.1 * db)


def pow2db(p):
    """Power ratio to decibels"""
    return 10 * np.log10(p)


def enob(sinad):
    """Calculate ENOB from SINAD"""
    return (sinad - 1.76) / 6.02


def snr_theor(n):
    """Theoretical SNR of an ideal n-bit ADC in dB"""
    return 6.02 * n + 1.76


def noise_floor(snr, m):
    """Noise floor of the m-point FFT in dB"""
    return -snr - 10 * np.log10(m / 2)


def harmonics(psp, fft_n, ref_pow, sample_freq, leak=20, n=5, window='hanning'):
    """Obtain first n harmonics properties from power spectrum"""
    # Coherence Gain and Noise Power Bandwidth for different windows
    win_params = {'uniform': {'cg': 1.0, 'npb': 1.0},
                  'hanning': {'cg': 0.5, 'npb': 1.5},
                  'hamming': {'cg': 0.54, 'npb': 1.36},
                  'blackman': {'cg': 0.42, 'npb': 1.73}}[window]
    fft_n = len(psp) * 2  # one side spectrum provided
    df = sample_freq / fft_n
    # calculate fundamental frequency
    fund_bin = np.argmax(psp)
    fund_freq = (np.sum([psp[i] * i * df for i in range(fund_bin - leak, fund_bin + leak + 1)]) /
                 np.sum(psp[fund_bin - leak: fund_bin + leak + 1]))
    # calculate harmonics info
    h = []
    for i in range(1, n + 1):
        h_i = {'num': i}
        zone_freq = (fund_freq * i) % sample_freq
        h_i['freq'] = sample_freq - zone_freq if zone_freq >= (sample_freq / 2) else zone_freq
        h_i['central_bin'] = int(h_i['freq'] / df)
        h_i['bins'] = np.array(range(h_i['central_bin'] - leak, h_i['central_bin'] + leak + 1))
        h_i['pow'] = ((1 / win_params['cg']) ** 2) * np.sum(psp[h_i['bins']]) / win_params['npb']
        h_i['vrms'] = np.sqrt(h_i['pow'])
        if i == 1:
            h_i['db'] = '%.2f dBFS' % pow2db(h_i['pow'] / ref_pow)
        else:
            h_i['db'] = '%.2f dBc' % pow2db(h_i['pow'] / h[0]['pow'])
        h += [h_i]
    return h


def signal_noise(psp, harmonics):
    """Obtain different signal+noise characteristics from spectrum"""
    # noise + distortion power
    nd_psp = np.copy(psp)
    nd_psp[harmonics[0]['bins']] = 0  # remove main harmonic
    nd_psp[0] = 0  # remove dc
    nd_pow = sum(nd_psp)
    # noise power
    n_psp = np.copy(psp)
    for h in harmonics:
        n_psp[h['bins']] = 0  # remove all harmonics
    n_psp[0] = 0  # remove dc
    n_pow = sum(n_psp)
    # distortion power
    d_pow = np.sum([h['pow'] for h in harmonics]) - harmonics[0]['pow']
    # calculate results
    sinad = pow2db(harmonics[0]['pow'] / nd_pow)
    thd = pow2db(harmonics[0]['pow'] / d_pow)
    snr = pow2db(harmonics[0]['pow'] / n_pow)
    sfdr = pow2db(max(nd_psp) / harmonics[0]['pow'])
    return sinad, thd, snr, sfdr


def analyze(sig, adc_bits, adc_vref, adc_freq, window='hanning'):
    """Do spectral analysis for ADC samples"""
    # Calculate some useful parameters
    sig_vpeak_max = adc_vref / 2
    sig_vrms_max = sig_vpeak_max / np.sqrt(2)
    sig_pow_max = sig_vrms_max ** 2
    ref_pow = sig_pow_max
    adc_prd = 1 / adc_freq
    adc_quants = 2 ** adc_bits
    dv = adc_vref / adc_quants
    sig_n = len(sig)
    dt = 1 / adc_freq
    fft_n = sig_n
    df = adc_freq / fft_n
    win_coef = {'uniform': np.ones(sig_n),
                'hanning': np.hanning(sig_n)}[window]
    sp_leak = 20  # spectru leak bins
    h_n = 5  # harmonics number

    # Convert samples to voltage
    sig_v = sig * dv

    # Remove DC and apply window
    sig_dc = np.mean(sig_v)
    sig_windowed = (sig_v - sig_dc) * win_coef

    # Calculate one-side amplitude spectrum (Vrms)
    asp = np.sqrt(2) * np.abs(np.fft.rfft(sig_windowed)) / sig_n

    # Calculate one-side power spectrum (Vrms^2)
    psp = np.power(asp, 2)
    psp_db = pow2db(psp / ref_pow)

    # Calculate harmonics
    h = harmonics(psp=psp, fft_n=fft_n, ref_pow=ref_pow, sample_freq=adc_freq, leak=sp_leak, n=h_n, window=window)

    # Input signal parameters (based on 1st harmonic)
    sig_pow = h[0]['pow']
    sig_vrms = h[0]['vrms']
    sig_vpeak = sig_vrms * np.sqrt(2)
    sig_freq = h[0]['freq']
    sig_prd = 1 / sig_freq

    # Calculate SINAD, THD, SNR, SFDR
    adc_sinad, adc_thd, adc_snr, adc_sfdr = signal_noise(psp, h)

    # Calculate ENOB
    # sinad correction to normalize ENOB to full-scale regardless of input signal amplitude
    adc_enob = enob(adc_sinad + pow2db(ref_pow / sig_pow))

    # Calculate Noise Floor
    adc_noise_floor = noise_floor(adc_snr, fft_n)

    # Create plots
    fig = plt.figure(figsize=(14, 7))
    gs = matplotlib.gridspec.GridSpec(2, 2, width_ratios=[3, 1])

    # Time plot
    ax_time = plt.subplot(gs[0, 0])
    ax_time_xlim = min(sig_n, int(5 * sig_prd / dt))
    ax_time.plot(np.arange(0, ax_time_xlim), sig[:ax_time_xlim], color='C0')
    ax_time.set(ylabel='ADC code', ylim=[0, adc_quants])
    ax_time.set(yticks=list(range(0, adc_quants, adc_quants // 8)) + [adc_quants - 1])
    ax_time.set(xlabel='Sample', xlim=[0, ax_time_xlim - 1])
    ax_time.set(xticks=range(0, ax_time_xlim, max(1, ax_time_xlim // 20)))
    ax_time.grid(True)
    ax_time_xsec = ax_time.twiny()
    ax_time_xsec.set(xticks=ax_time.get_xticks())
    ax_time_xsec.set(xbound=ax_time.get_xbound())
    ax_time_xsec.set_xticklabels(['%.02f' % (x * dt * 1e3) for x in ax_time.get_xticks()])
    ax_time_xsec.set_xlabel('Time, ms')
    ax_time_ysec = ax_time.twinx()
    ax_time_ysec.set(yticks=ax_time.get_yticks())
    ax_time_ysec.set(ybound=ax_time.get_ybound())
    ax_time_ysec.set_yticklabels(['%.02f' % (x * dv) for x in ax_time.get_yticks()])
    ax_time_ysec.set_ylabel('Voltage, V')

    # Frequency plot
    ax_freq = plt.subplot(gs[1, 0])
    ax_freq.plot(np.arange(0, len(psp_db)), psp_db, color='C0', zorder=0, label="Spectrum")
    for h_i in h:
        ax_freq.text(h_i['central_bin'] + 2, psp_db[h_i['central_bin']], str(h_i['num']),
                     va='bottom', ha='left', weight='bold')
        ax_freq.plot(h_i['bins'], psp_db[h_i['bins']], color='C4')
    ax_freq.plot(0, 0, color='C4', label="Harmonics")
    ax_freq.set(ylabel='dB', ylim=[-150, 10])
    ax_freq.set(xlabel='Sample', xlim=[0, fft_n / 2])
    ax_freq.set(xticks=list(range(0, fft_n // 2, fft_n // 32)) + [fft_n // 2 - 1])
    ax_freq.grid(True)
    ax_freq.legend(loc="lower right", ncol=3)
    ax_freq_sec = ax_freq.twiny()
    ax_freq_sec.set_xticks(ax_freq.get_xticks())
    ax_freq_sec.set_xbound(ax_freq.get_xbound())
    ax_freq_sec.set_xticklabels(['%.02f' % (x * df * 1e-3) for x in ax_freq.get_xticks()])
    ax_freq_sec.set_xlabel('Frequency, kHz')

    # Information plot
    ax_info = plt.subplot(gs[:, 1])
    ax_info.set(xlim=[0, 10], xticks=[], ylim=[0, 10], yticks=[])
    harmonics_str = '\n'.join(['%d%s @ %-10s : %s' % (h_i['num'], ['st', 'nd', 'rd', 'th', 'th'][h_i['num'] - 1],
                                                      '%0.3f kHz' % (h_i['freq'] * 1e-3),
                                                      h_i['db']) for h_i in h])
    ax_info_str = """
========= FFT ==========
Points           : {fft_n}
Freq. resolution : {fft_res:.4} Hz
Window           : {fft_window}

======= Harmonics ======
{harmonics_str}

===== Input signal =====
Frequency        : {sig_freq:.4} kHz
Amplitude (Vpeak): {sig_vpeak:.4} V
DC offset        : {sig_dc:.4} V

========= ADC ==========
Sampling freq.   : {adc_freq:.4} kHz
Sampling period  : {adc_prd:.4} us
Reference volt.  : {adc_vref:.4} V
Bits             : {adc_bits} bits
Quants           : {adc_quants}
Quant            : {adc_quant:.4} mV
SNR              : {adc_snr:.4} dB
SINAD            : {adc_sinad:.4} dB
THD              : {adc_thd:.4} dB
ENOB             : {adc_enob:.4} bits
SFDR             : {adc_sfdr:.4} dBc
Noise floor      : {adc_nfloor:.4} dBFS
""".format(fft_n=fft_n,
           fft_res=df,
           fft_window=window,
           harmonics_str=harmonics_str,
           sig_freq=sig_freq * 1e-3,
           sig_vpeak=sig_vpeak,
           sig_dc=sig_dc,
           adc_freq=adc_freq * 1e-3,
           adc_prd=adc_prd * 1e6,
           adc_vref=adc_vref,
           adc_bits=adc_bits,
           adc_quants=adc_quants,
           adc_quant=dv * 1e3,
           adc_snr=adc_snr,
           adc_thd=adc_thd,
           adc_sinad=adc_sinad,
           adc_enob=adc_enob,
           adc_sfdr=adc_sfdr,
           adc_nfloor=adc_noise_floor)
    ax_info.text(1, 9.5, ax_info_str, va='top', ha='left', family='monospace')

    # General plotting settings
    plt.tight_layout()
    plt.style.use('bmh')

    # Show the result
    plt.show()
