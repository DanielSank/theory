#!/usr/bin/python

filenames = [
    'aliasing_real',
    'aliasing',
    'aliasing_real_spectra',
    'aliasing_sampling',
    'bands'
]


command = ''
for name in filenames:
    command += "echo {}.svg --export-pdf={}.pdf;\n".format(name, name)
command = "({}) |\nDISPLAY= inkscape --shell".format(command)


import os
os.system(command)

