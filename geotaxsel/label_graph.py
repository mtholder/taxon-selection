#!/usr/bin/env python
import sys
from .logs import info


def int_to_bit_set(int_coded):
    place = 1
    ret = set()
    while int_coded >= place:
        if (place & int_coded) != 0:
            ret.add(place)
        place *= 2
    return ret


class CCResolution(object):
    def __init__(self, subsets, sum_score):
        self.subsets = subsets
        self._score = sum_score

    @property
    def num_subsets(self):
        return len(self.subsets)

    @property
    def score(self):
        return self._score

    @property
    def subsets_as_bit_sets(self):
        bs_list = []
        for int_el in self.subsets:
            bs_list.append(int_to_bit_set(int_el))
        return bs_list


def cc_for_subset(leaves, possible_sets, par_cc):
    assert possible_sets
    assert par_cc is not None
    psl = 0
    for po in possible_sets:
        psl |= po
        if psl == leaves:
            break
    if psl != leaves:
        return None
    cache = par_cc.cache
    cache_hit = cache.get(leaves)
    if cache_hit is not None:
        # info(f"Cache hit for {len(leaves)} leaves and {len(possible_sets)} sets!")
        return cache_hit
    nsw = {}
    for ps in possible_sets:
        nsw[ps] = par_cc.b_subset_wts[ps]
    cc_cc = par_cc.compat_cache
    sub_cc = ConnectedComponent(b_leaves=leaves, b_subset_wts=nsw, top_cc=par_cc.top_cc)
    sub_cc.indent = "  " + par_cc.indent
    sub_cc.fill_resolutions()
    cache[leaves] = sub_cc
    return sub_cc


class ConnectedComponent(object):
    def __init__(
        self,
        label_set=None,
        freq=None,
        leaves=None,
        subset_wts=None,
        top_cc=None,
        b_leaves=None,
        b_subset_wts=None,
    ):
        self.str_leaves = set()
        self.resolutions = None
        self.min_num_subs = None
        self.max_num_subs = None
        self.indent = ""
        if top_cc is None:
            self.top_cc = self
            self.cache = {}
            self.compat_cache = {}
        else:
            self.top_cc = top_cc
            self.cache = top_cc.cache
            self.compat_cache = top_cc.compat_cache

        self.b_leaves = None
        self.b_subset_wts = None
        if leaves is not None:
            assert subset_wts is not None
            assert label_set is None
            assert freq is None
            assert b_leaves is None
            assert b_subset_wts is None
            self.str_subset_wts = dict(subset_wts)
            self.str_leaves.update(leaves)
        elif label_set is not None:
            assert subset_wts is None
            assert freq is not None
            assert b_leaves is None
            assert b_subset_wts is None
            self.str_subset_wts = {label_set: freq}
            self.str_leaves.update(label_set)
        else:
            assert subset_wts is None
            assert freq is None
            assert b_leaves is not None
            assert b_subset_wts is not None
            self.b_leaves = b_leaves
            self.b_subset_wts = b_subset_wts
        self.str_to_bit = None
        self.b_to_str = None

    def add_set(self, label_set, freq):
        self.str_leaves.update(label_set)
        self.str_subset_wts[label_set] = freq

    def write(self, out):
        if not self.resolutions:
            self.fill_resolutions()
        out.write(f"{len(self.str_leaves)} labels, {len(self.b_subset_wts)} subsets.\n")
        for x in range(self.min_num_subs, 1 + self.max_num_subs):
            res = self.resolutions.get(x)
            if res is not None:
                ts = self.to_str_subset_list(res.subsets_as_bit_sets)
                out.write(f"  {x} subsets, best score={res.score}: {ts}\n")
            else:
                out.write(f"  {x} subsets: no resolution\n")

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return self.leaves == other.leaves and self.subset_wts == other.subset_wts

    def to_str_subset_list(self, b_subset_list):
        ret = []
        for bit_set in b_subset_list:
            str_set = set()
            for el in bit_set:
                str_set.add(self.b_to_str[el])
            ret.append(str_set)
        return ret

    def to_int_encoding(self, str_subset):
        int_encoding = 0
        for label in self.str_leaves:
            int_encoding += self.str_to_bit[label]
        return int_encoding

    def add_resolution(self, res):
        sc = res.score
        num_sub = res.num_subsets
        if self.min_num_subs is None or num_sub < self.min_num_subs:
            self.min_num_subs = num_sub
        if self.max_num_subs is None or num_sub > self.max_num_subs:
            self.max_num_subs = num_sub
        prev = self.resolutions.get(num_sub)
        if prev is None or sc > prev.score:
            self.resolutions[num_sub] = res
            return True
        return False

    def _create_str_to_bit_encoding(self):
        self.b_leaves = 0
        self.b_subset_wts = {}
        self.str_to_bit = {}
        self.b_to_str = {}
        curr_place = 1
        for label in self.str_leaves:
            self.str_to_bit[label] = curr_place
            self.b_to_str[curr_place] = label
            self.b_leaves += curr_place
            curr_place *= 2
        for label_set, weight in self.str_subset_wts.items():
            self.b_subset_wts[self.to_int_encoding(label_set)] = weight

    def fill_resolutions(self):
        if self.cache is None:
            self.cache = {}
            assert self.compat_cache is None
            self.compat_cache = {}
        self.resolutions = {}
        if self.b_leaves is None:
            self._create_str_to_bit_encoding()
        else:
            assert self.b_subset_wts is not None
        all_subsets = list(self.b_subset_wts.keys())
        assert len(all_subsets) > 0
        one_subset = all_subsets[0]
        if len(all_subsets) == 1:
            assert one_subset == self.b_leaves
            sum_sc = self.b_subset_wts[one_subset]
            res = CCResolution(subsets=[one_subset], sum_score=sum_sc)
            self.add_resolution(res)
            return
        # forced, force_inc_leaves, sum_sc = self._force_maj_rule()
        # if force_inc_leaves == self.leaves:
        #     res = CCResolution(subsets=forced, sum_score=sum_sc)
        #     self.add_resolution(res)
        #     return
        # possible = [i for i in all_subsets if force_inc_leaves.isdisjoint(i)]
        forced = []
        force_inc_bits = 0
        sum_sc = 0
        possible = all_subsets  # alias due to trying out forcing maj-rule
        try:
            self._rec_fill_res(forced, force_inc_bits, possible, sum_sc)
        except AssertionError:
            print(f"fill_resolutions_possible = {possible}")
            print(f"fill_resolutions_forced = {forced}")
            raise

    def _rec_fill_res(self, in_subsets, in_labels, possible, sum_sc):
        one_poss_label = self._choose_one_possible_label(in_labels, possible)
        alternatives, others = [], set()
        for subset in possible:
            if one_poss_label & subset:
                alternatives.append(subset)
            else:
                others.add(subset)

        for idx, i in enumerate(alternatives):
            leaves_held = in_labels
            leaves_held |= i  # add label "i"
            i_sc = self.b_subset_wts[i]
            leaves_needed = self.b_leaves - leaves_held
            if not leaves_needed:
                ires = CCResolution(in_subsets + [i], i_sc + sum_sc)
                self.add_resolution(ires)
                continue
            poss_other = others.intersection(self.subsets_compat_with(i))
            # info(
            #     f"{self.indent}idx={idx} of #alt={len(alternatives)} |leaves|={len(self.leaves)} #others={len(others)}, #po_list={len(po_list)}"
            # )
            if not poss_other:
                continue
            sub_cc = cc_for_subset(leaves_needed, poss_other, self)
            if sub_cc is None:
                continue
            for num_subs, res in sub_cc.resolutions.items():
                bigger_res = CCResolution(
                    in_subsets + [i] + res.subsets, i_sc + sum_sc + res.score
                )
                self.add_resolution(bigger_res)

    def _choose_one_possible_label(self, in_leaves, possible):
        assert len(possible) > 0
        labels_left = self.b_leaves - in_leaves
        USE_MOST_FREQ = False
        if USE_MOST_FREQ:
            most_freq_label, most_freq_count = None, 0
            for label in int_to_bit_set(labels_left):
                reps = sum([1 for i in possible if label & i])
                if reps > most_freq_count:
                    most_freq_label, most_freq_count = label, reps
            return most_freq_label
        least_freq_label, least_freq_count = None, 1 + len(possible)
        for label in int_to_bit_set(labels_left):
            reps = sum([1 for i in possible if label & i])
            if reps < least_freq_count:
                least_freq_label, least_freq_count = label, reps
        return least_freq_label

    def subsets_compat_with(self, el):
        cs = self.compat_cache.get(el)
        # cs = None
        if cs is None:
            cs = frozenset([o for o in self.top_cc.b_subset_wts.keys() if not el & o])
            self.compat_cache[el] = cs
        return cs


class LabelGraph(object):
    def __init__(self):
        self.full_label_set = set()
        self.by_freq = []
        self.leaf_to_comp = {}
        self.components = set()

    def add_set(self, ls, freq):
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
            first_comp = ConnectedComponent(label_set=ls, freq=freq)
            self.components.add(first_comp)
            # debug(f"Adding comp {first_comp}")
        for el in ls:
            self.leaf_to_comp[el] = first_comp

    def merge_comps(self, first_comp, other):
        first_comp.leaves.update(other.leaves)
        first_comp.subset_wts.update(other.subset_wts)
        # debug(f"Removing comp {other} after merge with {first_comp}")
        self.components.remove(other)
        for el in other.leaves:
            self.leaf_to_comp[el] = first_comp

    def _get_sortable_comp_info(self):
        sortable = [
            (len(i.leaves), len(i.subset_wts), i.leaves, i) for i in self.components
        ]
        sortable.sort()
        return sortable

    def write(self, out):
        out.write(f"{len(sortable)} components:\n")
        for ind, el in enumerate(self._get_sortable_comp_info()):
            out.write(f"Component #{1 + ind}: ")
            el[-1].write(out)

    def write_components(self, fprefix):
        for ind, el in enumerate(self._get_sortable_comp_info()):
            fp = f"{fprefix}-{1+ind}.py"
            comp = el[-1]
            with open(fp, "w") as outp:
                info(f"Serializing component #{1 + ind} to {fp} ...")
                outp.write(
                    f"""#!/bin/env python
from geotaxsel.label_graph import ConnectedComponent
import sys

cc_{1+ind} = ConnectedComponent(leaves={repr(comp.leaves)},
                               subset_wts={repr(comp.subset_wts)}
                               )
cc_{1+ind}.write(sys.stdout)
"""
                )


def choose_most_common(rep_selections, num_to_select, num_trees):
    with open("cruft/TEMP_REP_SELS.python", "w") as outp:
        outp.write("x = ")
        outp.write(repr(rep_selections))
        outp.write("\n")
    lg = LabelGraph()

    for k, v in rep_selections.items():
        lg.add_set(k, v)
    lg.write_components("cruft/comp")
    sys.exit("early in choose_most_common\n")
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
