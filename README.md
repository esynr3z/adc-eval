# ADC evaluation tool

Tiny tools collection (Python [NumPy](https://numpy.org/)+[Matplotlib](https://matplotlib.org/) based) to do spectral analysis and calculate the key performance parameters of an ADC. Just collect some data from the ADC, specify basic ADC parameters and run analysis. See [example.ipynb](example.ipynb) (you will need [Jupyter Notebook](https://jupyter.org/) to be installed).

![analyser](analyser.png)

## CSV Evaluation Script

The included script `eval_csv.py` can be used to evaluate ADC samples stored in a CSV where each column is one "channel" or set of samples to be analized with one conversion per row.

Example:

|  set1  |  set2  |
|--------|--------|
| convr1 | convr1 |
| convr2 | convr2 |
| convr3 | convr3 |
| convr4 | convr4 |

## References
- [Analog Devices MT-003 TUTORIAL "Understand SINAD, ENOB, SNR, THD, THD + N, and SFDR so You Don't Get Lost in the Noise Floor"](https://www.analog.com/media/en/training-seminars/tutorials/MT-003.pdf)
- [National Instruments Application Note 041 "The Fundamentals of FFT-Based Signal Analysis and Measurement"](http://www.sjsu.edu/people/burford.furman/docs/me120/FFT_tutorial_NI.pdf)

Inspired by Linear Technology (now Analog Devices) [PScope](https://www.analog.com/en/technical-articles/pscope-basics.html) tool.

![pscope](pscope.png)
Image source: [Creating an ADC Using FPGA Resources WP - Lattice](https://www.latticesemi.com/-/media/LatticeSemi/Documents/WhitePapers/AG/CreatingAnADCUsingFPGAResources.ashx?document_id=36525)
