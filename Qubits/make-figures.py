#!/usr/bin/python

filenames = [
    'singleCircuitWithDrive',
    'coupledCircuits_capacitive',
    'coupledCircuits_inductive'
]


command = ''
for name in filenames:
    command += "echo {}.svg --export-pdf={}.pdf;\n".format(name, name)
command = "({}) |\nDISPLAY= inkscape --shell".format(command)


import os
os.system(command)

