#!/usr/bin/env pyth
import sys
from geotaxsel import read_taxonomy_stream


def main(taxonomy_fp=None):
    if taxonomy_fp is None:
        read_taxonomy_stream(sys.stdin)
    else:
        with open(taxonomy_fp, "r") as inp:
            read_taxonomy_stream(inp)


if __name__ == "__main__":
    main(None if len(sys.argv) == 1 else sys.argv[1])
