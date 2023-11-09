#!/usr/bin/env python
import sys
from tempfile import mkdtemp
import os
import subprocess
from .logs import info
from .label_graph import LabelGraph
import json


class ResolutionWrapper(object):
    def __init__(self, res_list):
        assert isinstance(res_list, list)
        self.min_num = None
        self.max_num = None
        self.by_size = {}
        for res in res_list:
            score = float(res["score"])
            assert isinstance(score, float)
            size = res["size"]
            assert isinstance(size, int)
            subsets_list = res["subsets"]
            fz_list = []
            for i in subsets_list:
                assert isinstance(i, list)
                fz_list.append(frozenset(i))
            if self.min_num is None or size < self.min_num:
                self.min_num = size
            if self.max_num is None or size > self.max_num:
                self.max_num = size
            self.by_size[size] = (score, fz_list)

    @property
    def size_width(self):
        return 1 + self.max_num - self.min_num


PROB_FN = "problems.csv"


def serialize_problems_for_most_common_choice(rep_selections, temp_dir=None):
    if temp_dir is None:
        temp_dir = mkdtemp(prefix="taxsel-scratch-", dir=os.curdir)

    with open(os.path.join(temp_dir, "TEMP_REP_SELS.py"), "w") as outp:
        outp.write("x = ")
        outp.write(repr(rep_selections))
        outp.write("\n")
    lg = LabelGraph()
    for k, v in rep_selections.items():
        lg.add_set(k, v)
    pref = os.path.join(temp_dir, "comp")
    written = lg.write_components(pref)
    tmp_loc = os.path.join(temp_dir, f".{PROB_FN}")
    with open(tmp_loc, "w") as flagf:
        for line in written:
            flagf.write(f"{line}\n")
    final_loc = os.path.join(temp_dir, PROB_FN)
    os.rename(tmp_loc, final_loc)
    return temp_dir


def _run_solver(inp_fp, out_fp, max_secs_per_run, num_greedy=0):
    hide_out = out_fp + ".HIDE"
    err_fp = out_fp + "-err.txt"
    invoc = ["max-weight-partition", inp_fp]
    if num_greedy > 0:
        invoc.append(str(num_greedy))
    done = False
    info(f"  Running: {' '.join(invoc)}")
    with open(err_fp, "w") as err_fo:
        with open(hide_out, "w") as out_fo:
            proc = subprocess.Popen(invoc, stdout=out_fo, stderr=err_fo)
            try:
                rc = proc.wait(timeout=max_secs_per_run)
            except subprocess.TimeoutExpired:
                proc.kill()
            else:
                if rc != 0:
                    raise RuntimeError(f"Invocation of {invoc} failed")
                done = True
    if done:
        os.rename(hide_out, out_fp)
    else:
        info(
            f"Solver did not complete on {inp_fp} within {max_secs_per_run}, trying with num_greedy steps set to {1 + num_greedy}"
        )
        _run_solver(inp_fp, out_fp, max_secs_per_run, num_greedy=1 + num_greedy)


def _run_solver_on_all(to_do_list, max_secs_per_run):
    for inp, outp in to_do_list:
        _run_solver(inp, outp, max_secs_per_run=max_secs_per_run)


def _ensure_problems_solved(inp_files, max_secs_per_run):
    to_do_list = []
    all_out_files = []
    for fs in inp_files:
        assert fs.endswith(".csv")
        stem = fs[:-4]
        expected_out = f"{stem}.json"
        all_out_files.append(expected_out)
        if not os.path.isfile(expected_out):
            to_do_list.append((fs, expected_out))

    if to_do_list:
        _run_solver_on_all(to_do_list, max_secs_per_run=max_secs_per_run)
    return all_out_files


def _process_resolution_files(resolution_files):
    obj_list = []
    for fp in resolution_files:
        with open(fp, "r") as inp:
            jobj = json.load(inp)
        obj_list.append(ResolutionWrapper(jobj))
    return obj_list


def choose_most_common(num_to_select, scratch_dir, max_secs_per_run=6000):
    prob_list_fp = os.path.join(scratch_dir, PROB_FN)

    with open(prob_list_fp, "r") as inp:
        inp_files = [i.strip() for i in inp]

    resolution_files = _ensure_problems_solved(
        inp_files, max_secs_per_run=max_secs_per_run
    )

    res_wrap_list = _process_resolution_files(resolution_files)
    # Sort by the ones with the smallest variation in size first
    sortable = [(i.size_width, i.min_num, id(i), i) for i in res_wrap_list]
    sortable.sort()
    in_order = [i[-1] for i in sortable]
    total_min = 0
    total_max = 0
    for res in in_order:
        total_max += res.max_num
        total_min += res.min_num
        print(res.min_num, res.max_num, total_min, total_max)

    if num_to_select < total_min or num_to_select > total_max:
        raise RuntimeError(
            f"Cannot select {num_to_select} lineages solutions range in [{total_min}, {total_max}]"
        )

    raise NotImplementedError("Not implemented yet")
    # with open("cruft/TEMP_REP_SELS.python", "w") as outp:
    #     outp.write("x = ")
    #     outp.write(repr(rep_selections))
    #     outp.write("\n")
    # lg = LabelGraph()
    #
    # for k, v in rep_selections.items():
    #     lg.add_set(k, v)
    # lg.write_components("cruft/comp")
    # sys.exit("early in choose_most_common\n")
    # lg.write(sys.stdout)

    # crs = dict(rep_selections)
    # by_freq.sort(reverse=True)
    # sel_set = set()
    # sel_leaf_set = set()
    # for times_sel, label_set in by_freq:
    #     if not label_set.isdisjoint(sel_leaf_set):
    #         info(f"  skipping label_set seen {times_sel} times: {label_set}")
    #         continue
    #     sel_set.add(label_set)
    #     sel_leaf_set.update(label_set)
    #     if len(sel_set) == num_to_select:
    #         if sel_leaf_set == full_label_set:
    #             return sel_set
    #         missing = full_label_set.difference(sel_leaf_set)
    #         sys.exit(f"Got to {len(sel_set)} while missing {missing}\n")
    #     if sel_leaf_set == full_label_set:
    #         sys.exit(f"Covered all leaves in only {len(sel_set)} sets\n")
    # assert False
