#!/usr/bin/env pyth
import sys
from geotaxsel import read_taxonomy_stream


def main(taxonomy_fp=None):
    if taxonomy_fp is None:
        clades = read_taxonomy_stream(sys.stdin)
    else:
        with open(taxonomy_fp, "r") as inp:
            clades = read_taxonomy_stream(inp)
    out = sys.stdout
    for line in clades:
        out.write(f"{line[0]}\t{line[1]}\n")


if __name__ == "__main__":
    main(None if len(sys.argv) == 1 else sys.argv[1])
