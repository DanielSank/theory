#!/usr/bin/python

filenames = [
    'coupled_circuits',
    'coupled_lines',
    'drive',
    'avoided_crossing',
]


command = ''
for name in filenames:
    command += "echo {}.svg --export-pdf={}.pdf;\n".format(name, name)
command = "({}) |\nDISPLAY= inkscape --shell".format(command)


import os
os.system(command)

