import argparse

from pmd.postprocessing import calculate_Tg


def main(args=None):
    args = parse_command_line(args)

    # TODO: setup a nice dict with all property options and
    # corresponding functions
    if args.property == 'Tg':
        calculate_Tg(args.result)


def parse_command_line(args=None):
    parser = argparse.ArgumentParser(description=(
        'Analyze simulation result files using the Pmd Analysis module.'))
    parser.add_argument('result', help='Result file')
    parser.add_argument('-p', '--property', type=str, help='Property name')
    args = parser.parse_args(args=args)

    return args


if __name__ == '__main__':
    main()
