#! /usr/bin/env python3
import sys

VERBOSE = True


def set_verbose(v=True):
    global VERBOSE
    VERBOSE = v


def debug(m):
    if VERBOSE:
        sys.stderr.write(f"greedy_mmd debug: {m}\n")


def info(msg):
    sys.stderr.write(f"{msg}\n")
