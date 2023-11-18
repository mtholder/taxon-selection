#!/usr/bin/env python

import sys

to_count = {}

for fp in sys.argv[1:]:
    with open(fp, "r") as inp:
        for line in inp:
            ls = line.strip()
            to_count[ls] = 1 + to_count.get(ls, 0)


by_counts = []
for k, v in to_count.items():
    by_counts.append((v,k))

by_counts.sort(reverse=True)
for count, name in by_counts:
    print(name, count)
