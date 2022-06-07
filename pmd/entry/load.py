import argparse
from typing import Optional, Sequence

from pmd.core.Pmd import Pmd


def main(args: Optional[Sequence[str]] = None):
    args = parse_command_line(args)
    output_dir = '.'
    if args.output_dir:
        output_dir = args.output_dir
    Pmd.load_config(args.config, output_dir)


def parse_command_line(args: Optional[Sequence[str]] = None):
    parser = argparse.ArgumentParser(
        description='Create simulation files from Pmd config files.')
    parser.add_argument('config', help='configuration file')
    parser.add_argument('-o',
                        '--output_dir',
                        type=str,
                        help='Output directory')
    args = parser.parse_args(args=args)

    return args


if __name__ == '__main__':
    main()
