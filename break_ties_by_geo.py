#!/usr/bin/env python
import argparse
import sys
import os
import re


def main(cut_branches_fp, chosen_tax_fp, centroid_fp):
    pass


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Expecting 3 arguments: cut_branches_fp chosen_tax_fp centroid_fp")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
