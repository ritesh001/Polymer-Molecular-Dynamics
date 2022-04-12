import os


def build_dir(output_dir):
    try:
        os.mkdir(output_dir)
    except OSError:
        pass


def register_kwargs(self_kwargs, input_kwargs):
    for key in self_kwargs:
        if key in input_kwargs:
            self_kwargs[key] = input_kwargs.get(key)
