import argparse
from typing import Optional, Sequence

from pmd.postprocessing import calculate_diffusivity, calculate_Tg


def main(args: Optional[Sequence[str]] = None):
    args = parse_command_line(args)

    # TODO: setup a nice dict with all property options and
    # corresponding functions
    if args.property == 'Tg':
        calculate_Tg(args.result)
    elif args.property == 'D':
        if 'result' in args.keys():
            calculate_diffusivity(args.result)
        else:
            calculate_diffusivity()


def parse_command_line(args: Optional[Sequence[str]] = None):
    parser = argparse.ArgumentParser(description=(
        'Analyze simulation result files using the Pmd Analysis module.'))
    parser.add_argument('result', help='Result file or folder')
    parser.add_argument('-p', '--property', type=str, help='Property name')
    args = parser.parse_args(args=args)

    return args


if __name__ == '__main__':
    main()
