#!/bin/bash
for ((i=0; i<10; i++)) ; do
	grep -P '^1\t' cut_branches_${i}.csv | sed -E 's/1\t//' > singleton_${i}.csv
	tail -n600 out-${i}.txt >c${i}.txt
done
python count_repeats.py c*.txt > by_counts.txt
python count_repeats.py singleton_* > by_count_singleton.tsv
