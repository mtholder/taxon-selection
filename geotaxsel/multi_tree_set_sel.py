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

    def pre_cond_check(self):
        try:
            assert self.min_num == min(self.by_size.keys())
            assert self.max_num == max(self.by_size.keys())
        except:
            print(f"self.min_num = {self.min_num}")
            print(f"self.max_num = {self.max_num}")
            print(f"self.by_size = {self.by_size}")
            raise

    def absorb(self, other, final_n, min_num_after):
        # print(
        #     f" enter absorb [{self.min_num}, {self.max_num}], [{other.min_num}, {other.max_num}], {final_n}, {min_num_after}"
        # )
        self.pre_cond_check()
        other.pre_cond_check()
        nmin = self.min_num + other.min_num
        nmax = self.max_num + other.max_num
        cropped_nmax = min(nmax, final_n - min_num_after)

        o_min = other.min_num
        o_max = other.max_num
        nbs = {}
        for tot in range(nmin, 1 + cropped_nmax):
            best_score = float("-inf")
            best_pair_lists = None
            for this_idx in range(self.min_num, 1 + self.max_num):
                other_idx = tot - this_idx
                if other_idx < o_min:
                    break
                if other_idx > o_max:
                    continue
                s_el = self.by_size.get(this_idx)
                if s_el is None:
                    continue
                o_el = other.by_size.get(other_idx)
                if o_el is None:
                    continue
                s_score, s_subsets = s_el
                o_score, o_subsets = o_el
                curr_score = s_score + o_score
                if curr_score > best_score:
                    best_score = curr_score
                    best_pair_lists = (s_subsets, o_subsets)
            if best_pair_lists is not None:
                nbs[tot] = (best_score, best_pair_lists[0] + best_pair_lists[1])
            else:
                info(f" No size combo added to {tot} !!")
        self.by_size = nbs
        self.min_num = nmin
        self.max_num = cropped_nmax
        self.pre_cond_check()
        # print(f" exit absorb [{self.min_num}, {self.max_num}] by_size len={len(self.by_size)}")


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
        # print(res.min_num, res.max_num, total_min, total_max)

    if num_to_select < total_min or num_to_select > total_max:
        raise RuntimeError(
            f"Cannot select {num_to_select} lineages solutions range in [{total_min}, {total_max}]"
        )

    num_subproblems = len(in_order)
    min_size_after = [0] * num_subproblems
    cum_min = 0
    for idx, res in enumerate(reversed(in_order)):
        if idx > 0:
            neg_idx = -(1 + idx)
            min_size_after[neg_idx] = cum_min
        cum_min += res.min_num

    final = in_order[0]
    for idx, res in enumerate(in_order):
        if idx == 0:
            continue
        else:
            assert res is not final
        final.absorb(res, final_n=num_to_select, min_num_after=min_size_after[idx])
    return final.by_size[num_to_select]
