import argparse

import numpy as np
import spectrum
import csv

DEFAULT_ADC_FREQ = 20_000_000
DEFAULT_ADC_BITS = 12
DEFAULT_ADC_VREF = 3.3

def main():
    with open(args.input_file) as fp:
        data = tuple(map(tuple, csv.reader(fp)))

    for channel_data in zip(*data):
        channel_data = np.asarray(tuple(map(int, channel_data)))
        spectrum.analyze(
            channel_data,
            args.bits,
            args.vref,
            args.frequency,
            window='hanning'
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser('ADC S2N/ENOB Eval Tool')
    parser.add_argument(
        'input_file',
        help='The input file from which to draw samples.  Should be CSV with' \
            'each channel to analyze being its own column and each row a ' \
            'sample.'
    )
    parser.add_argument(
        '-f', '--frequency', default=DEFAULT_ADC_FREQ,
        help=f'The sampling frequency of the ADC. Defaults to {DEFAULT_ADC_FREQ}.'
    )
    parser.add_argument(
        '-v', '--vref', default=DEFAULT_ADC_VREF,
        help=f'The voltage ref for the ADC. Defaults to {DEFAULT_ADC_VREF}.'
    )
    parser.add_argument(
        '-b', '--bits', default=DEFAULT_ADC_BITS,
        help=f'The number of bits of resolution the ADC has. Defaults ot {DEFAULT_ADC_BITS}'
    )

    args = parser.parse_args()

    main()