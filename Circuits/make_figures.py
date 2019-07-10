#!/usr/bin/python

filenames = [
    'distributed_resonator_quarter_wave',
    'mutual_inductance',
    'parallel_series_equivalent',
    'loaded_mode'
]

import os

command = ''
for name in filenames:
    command += "echo {}/{}.svg --export-pdf={}/{}.pdf;\n".format(
            os.getcwd(), name, os.getcwd(), name)
command = "({}) |\nDISPLAY= inkscape --shell".format(command)


os.system(command)

