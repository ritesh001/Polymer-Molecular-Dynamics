import os
import sys
from typing import Callable, Tuple


def validate_options(cls, options: Tuple[str]):
    num_given_options = sum(cls.__dict__[f"_{opt}"] is not None
                            for opt in options)
    option_string = ", ".join(opt for opt in options)
    if num_given_options == 0:
        raise ValueError(f'One of {option_string} has to be provided '
                         f'for the {cls} object')
    elif num_given_options > 1:
        raise ValueError(f'Only one of {option_string} can be provided '
                         f'for the {cls} object')


def build_dir(func: Callable) -> Callable:

    def wrapper_build_dir(*args, **kwargs):
        if 'output_dir' in kwargs:
            os.makedirs(kwargs['output_dir'], exist_ok=True)
        elif len(args) > 1:
            os.makedirs(args[1], exist_ok=True)
        return func(*args, **kwargs)

    return wrapper_build_dir


class HiddenPrints:

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
