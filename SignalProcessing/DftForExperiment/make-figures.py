import os

filenames = [
        'baseband',
        'leakage',
        'power',
        ]

for filename in filenames:
    os.system(f"inkscape {filename}.svg --export-filename={filename}.pdf")
