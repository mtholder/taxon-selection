#!/usr/bin/env python
import sys
from .logs import info


class ConnectedComponent(object):
    def __init__(self, label_set, freq):
        self.leaves = set()
        self.leaves.update(label_set)
        self.subset_wts = {label_set: freq}

    def add_set(self, label_set, freq):
        self.leaves.update(label_set)
        self.subset_wts[label_set] = freq

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self.leaves == other.leaves and self.subset_wts == other.subset_wts


class LabelGraph(object):
    def __init__(self):
        self.full_label_set = set()
        self.by_freq = []
        self.leaf_to_comp = {}
        self.components = set()

    def add_set(self, ls, freq):
        info(f"adding one with {len(ls)} and freq={freq}")
        self.full_label_set.update(ls)
        self.by_freq.append((freq, ls))
        comp_set = set()
        for el in ls:
            component = self.leaf_to_comp.get(el)
            if component is not None:
                comp_set.add(component)
        if comp_set:
            comp_list = list(comp_set)
            first_comp = comp_list[0]
            for c in comp_list[1:]:
                self.merge_comps(first_comp, c)
            first_comp.add_set(ls, freq)
        else:
            first_comp = ConnectedComponent(ls, freq)
            self.components.add(first_comp)
            info(f"Adding comp {first_comp}")
        for el in ls:
            self.leaf_to_comp[el] = first_comp

    def merge_comps(self, first_comp, other):
        first_comp.leaves.update(other.leaves)
        first_comp.subset_wts.update(other.subset_wts)
        info(f"Removing comp {other} after merge with {first_comp}")
        self.components.remove(other)
        for el in other.leaves:
            self.leaf_to_comp[el] = first_comp

    def write(self, out):
        sortable = [
            (len(i.leaves), len(i.subset_wts), i.leaves) for i in self.components
        ]
        sortable.sort()
        out.write(f"{len(sortable)} components:\n")
        for el in sortable:
            out.write(f"  {el}\n")


def choose_most_common(rep_selections, num_to_select, num_trees):
    # with open("cruft/TEMP_REP_SELS.python", "w") as outp:
    #     outp.write("x = ")
    #     outp.write(repr(rep_selections))
    #     outp.write("\n")
    lg = LabelGraph()

    for k, v in rep_selections.items():
        lg.add_set(k, v)

    lg.write(sys.stdout)
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
