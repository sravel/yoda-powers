def build_parser():
    import argparse
    parser = argparse.ArgumentParser(description='Process some sdfsf.')
    parser.add_argument('integers', metavar='N', type=int, nargs='+',
                        help='An integer for the accumulator.')
    parser.add_argument('-i', '--identity', type=int, required=True, default=0,
                        help='the default result for no arguments '
                             '(default: 0)')
    parser.add_argument('--sum', dest='accumulate', action='store_const',
                        const=sum, default=max, required=True,
                        help='Sum the integers (default: find the max).')
    return parser


def main():
    args = build_parser().parse_args()
    print(args.accumulate(args.integers))


if __name__ == '__main__':
    main()
