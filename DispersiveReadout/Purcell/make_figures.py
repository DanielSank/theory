#!/usr/bin/python

filenames = [
    'screenshot_SPICE',
    'purcell_theory_vs_numerical'
]

import os

command = ''
for name in filenames:
    command += "echo {}/{}.svg --export-pdf={}/{}.pdf;\n".format(
            os.getcwd(), name, os.getcwd(), name)
command = "({}) |\nDISPLAY= inkscape --shell".format(command)


os.system(command)

