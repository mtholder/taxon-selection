#! /usr/bin/env python3
__version__ = "0.0.1a"  # sync with setup.py
from .logs import set_verbose, info, debug
from .taxonomy import CladeDef, Ranks, read_taxonomy_stream
from .geo_tree_parser import parse_geo_and_tree
from .greedy_mmd import ultrametric_greedy_mmd, greedy_mmd
