#!/usr/bin/env python
import sys


def stripped_strings(fp):
    names = []
    with open(fp, "r") as inp:
        names = [i.strip() for i in inp if i.strip()]
    return names


def main(singletons_fp, ranked_cuts_fp):
    singletons = stripped_strings(singletons_fp)
    ranked_strings = stripped_strings(ranked_cuts_fp)
    top_ranked = frozenset([i.split("\t")[0].strip() for i in ranked_strings])
    for s in singletons:
        if s not in top_ranked:
            sys.stderr.write(f"{s} not a top choice!\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Expecting 3 arguments: singletons_fp ranked_cuts_fp")
    main(sys.argv[1], sys.argv[2])
