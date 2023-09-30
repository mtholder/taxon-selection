#!/usr/bin/env python
import sys


def read_taxonomy_stream(inp):
    pass


def main(taxonomy_fp=None):
    if taxonomy_fp is None:
        read_taxonomy_stream(sys.stdin)
    else:
        with open(taxonomy_fp, "r") as inp:
            read_taxonomy_stream(inp)


if __name__ == "__main__":
    main(None if len(sys.argv) == 1 else sys.argv[1])
