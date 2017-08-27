#!/usr/bin/python

filenames = [
    'claw_coupler',
    'cosine_potential',
    'coupled_circuits',
    'LC_oscillator',
    'LC_oscillator_annotated',
    'normal_metal_states',
    'readout_resonator',
    'superconducting_metal_states',
    'tls',
    'transmon',
    'transmon_magnetic_noise',
    'wavefunctions',
]

command = ''
for name in filenames:
    command += "echo {}.svg --export-pdf={}.pdf;\n".format(name, name)
command = "({}) |\nDISPLAY= inkscape --shell".format(command)

import os
os.system(command)

