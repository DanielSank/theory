#!/usr/bin/python

filenames = [
    'figure1',
    'figure2',
    'figure3',
]

import os

command = ''
for name in filenames:
    command += "echo {}/{}.svg --export-pdf={}/{}.pdf;\n".format(
            os.getcwd(), name, os.getcwd(), name)
command = "({}) |\nDISPLAY= inkscape --shell".format(command)


os.system(command)

