#!/usr/bin/env python
import sys
from .logs import info


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


TEST_ANSWER = set(
    [
        frozenset({"Microgale pusilla"}),
        frozenset({"Limnogale mergulus"}),
        frozenset({"Microgale cowani"}),
        frozenset({"Microgale talazaci"}),
        frozenset({"Microgale gymnorhyncha"}),
        frozenset({"Microgale thomasi"}),
        frozenset({"Microgale drouhardi"}),
        frozenset({"Microgale jenkinsae"}),
        frozenset({"Microgale parvula"}),
        frozenset({"Microgale brevicaudata"}),
        frozenset({"Microgale grandidieri"}),
        frozenset({"Microgale dobsoni"}),
        frozenset({"Microgale taiva"}),
        frozenset({"Microgale fotsifotsy"}),
        frozenset({"Microgale dryas"}),
        frozenset({"Microgale longicaudata"}),
        frozenset({"Microgale nasoloi"}),
        frozenset(
            {
                "Microgale gracilis",
                "Microgale monticola",
                "Microgale majori",
                "Microgale principula",
                "Microgale jobihely",
            }
        ),
        frozenset({"Microgale soricoides"}),
    ]
)


def all_in_test(prev, ne):
    for el in prev:
        if el not in TEST_ANSWER:
            return False
    return ne in TEST_ANSWER


class ConnectedComponent(object):
    def __init__(self, label_set=None, freq=None, leaves=None, subset_wts=None):
        self.leaves = set()
        self.resolutions = {}
        self.min_num_subs = None
        self.max_num_subs = None

        if label_set is None:
            assert leaves is not None
            assert subset_wts is not None
            self.subset_wts = dict(subset_wts)
            self.leaves.update(leaves)
        else:
            assert freq is not None
            self.subset_wts = {label_set: freq}
            self.leaves.update(label_set)

    def add_set(self, label_set, freq):
        self.leaves.update(label_set)
        self.subset_wts[label_set] = freq

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

    def _get_one_label(self):
        one_subset = next(iter(self.subset_wts.keys()))
        return next(iter(one_subset))

    def _count_trees(self):
        num_trees = 0
        one_label = self._get_one_label()
        for subset, score in self.subset_wts.items():
            if one_label in subset:
                num_trees += score
        return num_trees

    def fill_resolutions(self):
        self.resolutions = {}
        all_subsets = list(self.subset_wts.keys())
        assert len(all_subsets) > 0
        one_subset = all_subsets[0]
        if len(all_subsets) == 1:
            assert one_subset == self.leaves
            res = CCResolution(
                subsets=[one_subset], sum_score=self.subset_wts[one_subset]
            )
            self.add_resolution(res)
            return
        force_inc_subsets = set()
        force_inc_leaves = set()
        sum_sc = 0

        FORCING_MAJ_RULE = False
        if FORCING_MAJ_RULE:
            num_trees = self._count_trees()
            info(f"num_trees = {num_trees}")
            # Force include any subsets in the majority
            maj_cutoff = 1 + (num_trees // 2)
            for subset, wt in self.subset_wts.items():
                if wt >= maj_cutoff:
                    force_inc_subsets.add(subset)
                    assert force_inc_leaves.isdisjoint(subset)
                    force_inc_leaves.update(subset)
                    sum_sc += wt
            info(
                f"Maj-rule forcing: {len(force_inc_subsets)} subsest, {len(force_inc_leaves)} labels"
            )
        forced = list(force_inc_subsets)
        if force_inc_leaves == self.leaves:
            res = CCResolution(subsets=forced, sum_score=sum_sc)
            self.add_resolution(res)
            return
        possible = [i for i in all_subsets if force_inc_leaves.isdisjoint(i)]
        # if len(possible) == 0:
        #     print("forced", force_inc_leaves)
        #     print("all", self.leaves)
        #     print("missing", self.leaves.difference(force_inc_leaves))
        #     for el in forced:
        #         print(f"  forced freq={self.subset_wts[el]} el=", el)
        #     for el in all_subsets:
        #         if el not in forced:
        #             print(f"  poss freq={self.subset_wts[el]} el=", el)
        try:
            self._rec_fill_res(forced, force_inc_leaves, possible, sum_sc)
        except AssertionError:
            print(f"top_possible = {possible}")
            print(f"top_forced = {forced}")
            raise

    def _rec_fill_res(self, in_sets, in_leaves, possible, sum_sc, prev_str=""):
        assert len(possible) > 0
        labels_left = self.leaves.difference(in_leaves)
        USE_MOST_FREQ = False
        if USE_MOST_FREQ:
            most_freq_label, most_freq_count = None, 0
            for label in labels_left:
                reps = sum([1 for i in possible if label in i])
                if reps > most_freq_count:
                    most_freq_label, most_freq_count = label, reps
            one_poss_label = most_freq_label
        else:
            least_freq_label, least_freq_count = None, 1 + len(possible)
            for label in labels_left:
                reps = sum([1 for i in possible if label in i])
                if reps < least_freq_count:
                    least_freq_label, least_freq_count = label, reps
            one_poss_label = least_freq_label
        alternatives, others = [], []
        for subset in possible:
            dest = alternatives if one_poss_label in subset else others
            dest.append(subset)
        # indent = "  " * (1 + len(in_sets))
        # info(
        #     f"{indent}_rfr({len(self.leaves) - len(in_leaves)} leaves more. {len(possible)} others)"
        # )

        for idx, i in enumerate(alternatives):
            next_str = f"{prev_str} ({idx} / {len(alternatives)})"
            info(next_str)
            po_list = [o for o in others if i.isdisjoint(o)]
            po_ls = set(i)
            assert po_ls.isdisjoint(in_leaves)
            complements_in_sets = False
            num_needed = len(self.leaves) - len(in_leaves)
            for o in po_list:
                po_ls.update(o)
                if len(po_ls) >= num_needed:
                    assert len(po_ls) == num_needed
                    complements_in_sets = True
                    break
            if not complements_in_sets:
                if all_in_test(in_sets, i):
                    info("  REJECTING branch due to failure to complement.")
                    print("in_sets =", in_sets)
                    print("i =", i)
                    print("others =", others)
                    print("po_list =", po_list)
                    assert False
                continue
            try:
                self._ind_rec_fill(
                    inc_leaves=set(in_leaves),
                    newest=i,
                    prev_subs=list(in_sets),
                    poss_others=po_list,
                    sum_sc=sum_sc,
                    prev_str=next_str,
                )
            except AssertionError:
                print("_rfr_in_sets =", in_sets)
                print("_rfr_i =", i)
                print("_rfr_others =", others)
                print("_rfr_po_list =", po_list)
                raise

    def _ind_rec_fill(
        self, inc_leaves, newest, prev_subs, poss_others, sum_sc, prev_str=""
    ):
        assert inc_leaves.isdisjoint(newest)
        inc_leaves.update(newest)
        sum_sc += self.subset_wts[newest]
        prev_subs.append(newest)
        if inc_leaves == self.leaves:
            return self.add_resolution(CCResolution(prev_subs, sum_sc))
        if not poss_others:
            if all_in_test(prev_subs, newest):
                assert False
            return False
        # indent = "  " * (1 + len(prev_subs))
        # info(
        #     f"{indent} _irfr({len(self.leaves) - len(inc_leaves)} leaves more. {len(poss_others)} others"
        # )
        one_other = poss_others[0]
        if len(poss_others) == 1:
            inc_leaves.update(one_other)
            if inc_leaves != self.leaves:
                if all_in_test(prev_subs, newest):
                    assert False
                return False
            prev_subs.append(one_other)
            sum_sc += self.subset_wts[one_other]
            return self.add_resolution(CCResolution(prev_subs, sum_sc))
        try:
            return self._rec_fill_res(
                prev_subs, inc_leaves, poss_others, sum_sc, prev_str=prev_str
            )
        except AssertionError:
            print("_irfr_prev_subs =", prev_subs)
            print("_irfr_poss_others =", poss_others)
            raise

        # others = [i for i in poss_others if inc_leaves.isdisjoint(i)]
        # ret = False
        # while others:
        #     next_new = others.pop(0)
        #     nr = self._continue_resolution(inc_leaves=set(inc_leaves),
        #                               newest=next_new,
        #                               prev_subs=list(prev_subs),
        #                               poss_others=others,
        #                               sum_sc=sum_sc)
        #     ret = ret or nr
        # return ret

    def write(self, out):
        if not self.resolutions:
            self.fill_resolutions()
        out.write(f"{len(self.leaves)} labels, {len(self.subset_wts)} subsets.\n")
        for x in range(self.min_num_subs, 1 + self.max_num_subs):
            res = self.resolutions.get(x)
            if res is not None:
                out.write(f"  {x} subsets, best score={res.score}: {res.subsets}\n")
            else:
                out.write(f"  {x} subsets: no resolution\n")

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
