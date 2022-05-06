import os, sys
from typing import Callable


def build_dir(func: Callable) -> Callable:

    def wrapper_build_dir(*args, **kwargs):
        output_dir_key = 'output_dir'
        output_dir = (kwargs[output_dir_key]
                      if output_dir_key in kwargs.keys() else args[1])
        os.makedirs(output_dir, exist_ok=True)
        return func(*args, **kwargs)

    return wrapper_build_dir


class HiddenPrints:

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
