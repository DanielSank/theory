#!/usr/bin/python

filenames = [
    'lumped_model',
]

import os

command = ''
for name in filenames:
    command += "echo {}/{}.svg --export-pdf={}/{}.pdf;\n".format(
            os.getcwd(), name, os.getcwd(), name)
command = "({}) |\nDISPLAY= inkscape --shell".format(command)


os.system(command)

