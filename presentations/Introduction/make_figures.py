#!/usr/bin/python

filenames = [
    'mech_OR',
    'NAND',
    'transistor',
    'transistor_superposition',
    'transistors_microscopic',
    'classical_vs_quantum_paths',
    'placeholder',
]

command = ''
for name in filenames:
    command += "echo {}.svg --export-pdf={}.pdf;\n".format(name, name)
command = "({}) |\nDISPLAY= inkscape --shell".format(command)

import os
os.system(command)

