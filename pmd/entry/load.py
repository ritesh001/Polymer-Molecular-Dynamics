import argparse
import pmd


def main(args=None):
    args = parse_command_line(args)
    # TODO: Make nicer import
    pmd.Pmd.load_config(args.config)


def parse_command_line(args=None):
    parser = argparse.ArgumentParser(
        description='Create simulation files from Pmd config files.')
    parser.add_argument('config', help='configuration file')
    args = parser.parse_args(args=args)

    return args


if __name__ == '__main__':
    main()
