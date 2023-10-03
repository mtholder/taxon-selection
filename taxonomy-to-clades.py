#!/usr/bin/env pyth
import sys
from geotaxsel import read_taxonomy_stream


def main(taxonomy_fp, out_cmw2mdd_fp):
    if taxonomy_fp is None:
        clades, cmw2curr = read_taxonomy_stream(sys.stdin)
    else:
        with open(taxonomy_fp, "r") as inp:
            clades, cmw2curr = read_taxonomy_stream(inp)
    out = sys.stdout
    for line in clades:
        out.write(f"{line[0]}\t{line[1]}\n")
    with open(out_cmw2mdd_fp, "w") as outp:
        outp.write("CMW_sciName\tsciName\n")
        sk = list(cmw2curr.keys())
        sk.sort()
        for k in sk:
            v = cmw2curr[k]
            outp.write(f"{k}\t{v}\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(
            "Expecting 2 arguments: (input) taxonomy.csv and (output) cmw-to-mdd-name-mapping.tsv"
        )
    main(sys.argv[1], sys.argv[2])
