#python3

import pathlib
import os


class change_dir:
    original: pathlib.Path
    target: pathlib.Path

    def __init__(self, directory: pathlib.Path) -> None:
        self.target = directory.absolute()
        self.original = pathlib.Path.cwd().absolute()

    def __enter__(self) -> None:
        os.chdir(self.target)

    def __exit__(self, type, value, traceback) -> None:
        del type, value, traceback
        os.chdir(self.original)


def main() -> None:
    for path in pathlib.Path(".").rglob("*.svg"):
        print(path)
        with change_dir(path.parent.absolute()):
            os.system(f"inkscape --export-filename={path.stem}.pdf {path.name}")


if __name__ == "__main__":
    main()
