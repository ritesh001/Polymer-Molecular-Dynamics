import os, sys


def build_dir(output_dir: str) -> None:
    try:
        os.mkdir(output_dir)
    except OSError:
        pass


class HiddenPrints:

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
