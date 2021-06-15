#!/usr/bin/python3

import os

filenames = [
    'aliasing_real',
    'aliasing',
    'aliasing_real_spectra',
    'aliasing_sampling',
    'bands'
]


for name in filenames:
    os.system(f'inkscape --export-type="pdf" {name}.svg')
