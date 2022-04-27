import os


def build_dir(output_dir: str) -> None:
    try:
        os.mkdir(output_dir)
    except OSError:
        pass
