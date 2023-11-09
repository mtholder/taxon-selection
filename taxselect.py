#!/usr/bin/env python
import argparse
import sys
import os
import re

from geotaxsel import (
    greedy_mmd,
    output_chosen_anc,
    parse_geo_and_tree,
    parse_geo,
    prune_taxa_without_sp_data,
    ultrametric_greedy_mmd,
    serialize_problems_for_most_common_choice,
    choose_most_common,
    PROB_FN,
    calc_dist,
)
import dendropy


class RunSettings(object):
    def __init__(
        self,
        country_name_fp=None,
        centroid_fp=None,
        tree_fp=None,
        name_mapping_fp=None,
        num_to_select=2,
        use_ultrametricity=True,
        clade_defs_fp=None,
        name_updating_fp=None,
        cut_branches_fp=None,
        tree_dir=None,
        ultrametric_tol=5e-5,
        scratch_dir=None,
        max_solver_seconds=6000,
    ):
        self.country_name_fp = country_name_fp
        self.centroid_fp = centroid_fp
        self.tree_fp = tree_fp
        self.name_mapping_fp = name_mapping_fp
        self.num_to_select = num_to_select
        self.use_ultrametricity = use_ultrametricity
        self.clade_defs_fp = clade_defs_fp
        self.name_updating_fp = name_updating_fp
        self.cut_branches_fp = cut_branches_fp
        self.tree_dir = tree_dir
        self.ultrametric_tol = ultrametric_tol
        self.scratch_dir = scratch_dir
        self.max_solver_seconds = max_solver_seconds


def record_clade_sel(sel, rep_selections):
    for anc in sel:
        leaves_below = list(anc.leaf_nodes())
        labels_below = frozenset([i.taxon.label for i in leaves_below])
        pn = rep_selections.get(labels_below, 0)
        rep_selections[labels_below] = 1 + pn


def create_most_common_groups_probs(
    geo_ret,
    centroid_fp=None,
    name_mapping_fp=None,
    num_to_select=1,
    use_ultrametricity=True,
    tree_dir=None,
    ultrametric_tol=5e-5,
):
    sp_by_name, clades, upham_to_iucn, new_names_for_leaves = geo_ret
    file_names = os.listdir(tree_dir)
    file_names.sort()
    sp_pat = re.compile(r"^([A-Z][a-z]+ +[-a-z0-9]+)$")
    rep_selections = {}
    for tree_n, el in enumerate(file_names):
        tree_fp = os.path.join(tree_dir, el)
        tree = dendropy.Tree.get(path=tree_fp, schema="newick")
        print(tree_fp)
        prune_taxa_without_sp_data(
            tree,
            frozenset(sp_by_name.keys()),
            upham_to_iucn=upham_to_iucn,
            name_mapping_fp=name_mapping_fp,
            centroid_fp=centroid_fp,
            clades=clades,
            new_names_for_leaves=new_names_for_leaves,
            sp_pat_in_tree=sp_pat,
        )
        if use_ultrametricity:
            sel = ultrametric_greedy_mmd(
                tree, num_to_select, sp_by_name, ultrametric_tol=ultrametric_tol
            )
        else:
            sel = greedy_mmd(tree, num_to_select, sp_by_name)
        record_clade_sel(sel, rep_selections)
    return rep_selections


def min_dist_to_member(loc, loc_set):
    min_d = float("inf")
    for i in loc_set:
        d = calc_dist(i, loc)
        if d < min_d:
            min_d = d
    return min_d


def choose_exemplars_by_geo_divergence(settings, geo_ret, final_subsets):
    chosen_labels = set()
    labels_to_check = set()
    label_to_group = {}
    num_to_select = len(final_subsets)
    for group in final_subsets:
        assert len(group) > 0
        if len(group) == 1:
            label = list(group)[0]
            chosen_labels.add(label)
        else:
            labels_to_check.update(group)
            for member in group:
                label_to_group[member] = group
    print(f"{len(chosen_labels)} forced choices due to singleton groups.")
    print(f"{len(labels_to_check)} other labels")
    # print(chosen_labels)
    if len(labels_to_check) == 0:
        return chosen_labels

    sp_by_name = geo_ret[0]
    # This version of the code only works when there is 1 location
    #       per taxon
    for k, v in sp_by_name.items():
        assert len(v.locations) == 1

    sp_2_loc = {}
    locs_chosen = set()
    for label, species in sp_by_name.items():
        loc = next(iter(species.locations))
        assert loc is not None
        if label in chosen_labels:
            locs_chosen.add(loc)
        else:
            sp_2_loc[label] = loc
    label_2_min_dist = {}
    max_min_d = float("-inf")
    next_chosen_label, next_chosen_loc = None, None
    for md_idx, pair in enumerate(sp_2_loc.items()):
        sp, loc = pair
        d = min_dist_to_member(loc, locs_chosen)
        label_2_min_dist[sp] = d
        if (md_idx + 1) % 100 == 0:
            sys.stderr.write(f" ... {sp} ==> min_dist {d}\n")
        if d > max_min_d:
            max_min_d = d
            next_chosen_label = sp
            next_chosen_loc = loc
    assert next_chosen_label is not None
    while True:
        print(
            f"Adding {next_chosen_label} with a min_dist of {max_min_d} from a previously chosen taxon"
        )
        chosen_labels.add(next_chosen_label)
        if len(chosen_labels) == num_to_select:
            break
        locs_chosen.add(next_chosen_loc)
        last_added_loc = next_chosen_loc
        group = label_to_group[next_chosen_label]
        for el in group:
            del label_2_min_dist[el]
        # update next_chosen...
        next_chosen_label, next_chosen_loc = None, None
        to_do = list(label_2_min_dist.keys())
        max_min_d = float("-inf")
        for sp in to_do:
            md = label_2_min_dist[sp]
            dist_to_most_recent = calc_dist(last_added_loc, sp_2_loc[sp])
            if dist_to_most_recent < md:
                label_2_min_dist[sp] = dist_to_most_recent
                md = dist_to_most_recent
            if md > max_min_d:
                max_min_d = md
                next_chosen_label = sp
                next_chosen_loc = loc
    group = label_to_group[next_chosen_label]
    assert len(group) == len(label_2_min_dist)

    for group in final_subsets:
        sg = set(group)
        iset = chosen_labels.intersection(sg)
        liset = len(iset)
        if liset != 1:
            print(f"group ({group}) has {liset} members in {chosen_labels}")
            assert liset == 1
    return chosen_labels


def run_tree_dir(settings):
    geo_ret = parse_geo(
        country_name_fp=settings.country_name_fp,
        centroid_fp=settings.centroid_fp,
        name_mapping_fp=settings.name_mapping_fp,
        clade_defs_fp=settings.clade_defs_fp,
        name_updating_fp=settings.name_updating_fp,
    )
    if settings.scratch_dir is not None:
        if not os.path.isdir(settings.scratch_dir):
            raise RuntimeError(f"scratch_dir '{settings.scratch_dir}' does not exist.")
        one_comp_py_fp = os.path.join(settings.scratch_dir, PROB_FN)
        need_most_common_prob = not os.path.isfile(one_comp_py_fp)
    else:
        need_most_common_prob = True
    if need_most_common_prob:
        rep_selections = create_most_common_groups_probs(
            geo_ret,
            centroid_fp=settings.centroid_fp,
            name_mapping_fp=settings.name_mapping_fp,
            num_to_select=settings.num_to_select,
            use_ultrametricity=settings.use_ultrametricity,
            tree_dir=settings.tree_dir,
            ultrametric_tol=settings.ultrametric_tol,
        )
        td = serialize_problems_for_most_common_choice(rep_selections)
    else:
        td = settings.scratch_dir
    final_sc, final_subsets = choose_most_common(
        num_to_select=settings.num_to_select,
        scratch_dir=td,
        max_secs_per_run=settings.max_solver_seconds,
    )
    output_chosen_anc(
        tree=None,
        cut_branches_fp=settings.cut_branches_fp,
        chosen_ancs=final_subsets,
    )
    taxa = choose_exemplars_by_geo_divergence(settings, geo_ret, final_subsets)
    sys.exit("early at the end of run_tree_dir")
    return 0


def run(settings):
    if settings.tree_dir is not None:
        return run_tree_dir(settings)
    sp_pat = re.compile("^([A-Z][a-z]+ +[-a-z0-9]+) [A-Z][A-Za-z]+ [A-Z]+$")
    assert tree_fp is not None
    tree, sp_by_name = parse_geo_and_tree(
        settings.country_name_fp,
        settings.centroid_fp,
        settings.name_mapping_fp,
        settings.tree_fp,
        clade_defs_fp=settings.clade_defs_fp,
        name_updating_fp=settings.name_updating_fp,
        sp_pat_in_tree=settings.sp_pat,
    )
    if settings.use_ultrametricity:
        sel = ultrametric_greedy_mmd(
            tree,
            settings.num_to_select,
            sp_by_name,
            ultrametric_tol=settings.ultrametric_tol,
        )
    else:
        sel = greedy_mmd(tree, settings.num_to_select, sp_by_name)
    output_chosen_anc(tree, settings.cut_branches_fp, sel)
    sys.exit("early exit\n")
    print("Selected:\n  {}\n".format("\n  ".join(sel)))
    return 0


def main():
    parser = argparse.ArgumentParser("taxselect.py")
    parser.add_argument(
        "--country-file",
        default=None,
        help="Optional list of countries, used if taxa have geo. "
        "listings as a sets of centroids within several countries",
    )
    parser.add_argument(
        "--name-mapping-file",
        default=None,
        help="Optional file tab-separated file with header and "
        '"name_in_tree" then "IUCN_Name" used in country-file as the '
        "first two columns of each line",
    )
    parser.add_argument(
        "--name-updating-file",
        default=None,
        help="Optional file tab-separated file with header and "
        '"CMW_sciName" then "sciName" where sciName may appear"'
        "in the clades file, but CMW_sciName may be used in the"
        "tree.",
    )
    parser.add_argument(
        "--centroid-file",
        default=None,
        required=True,
        help="Filepath to CSV file with centroid locations for species. "
        "If a country-file and name-mapping-file are being used, the "
        "program expects the fields to be: \n"
        "Species_no, binomial, Country_ID, Country, Longitude, Latitude\n"
        "If there is just one line per species, the field ordering should be:\n"
        "name,x,y\n",
    )
    parser.add_argument(
        "--clade-defs-file",
        default=None,
        required=False,
        help="Filepath to tab-separated file no header that is the output of "
        "taxonomy-to-clades.py. Each line should contain a name in the first "
        "column and clade definition in the second row.",
    )
    parser.add_argument(
        "--cut-branches-file",
        default=None,
        required=False,
        help="Optional filepath to output CSV describing the branch cut points.",
    )
    parser.add_argument(
        "--tree-file",
        default=None,
        required=False,
        help="path to NEXUS file with a single ultrametric tree",
    )
    parser.add_argument(
        "--tree-dir",
        default=None,
        required=False,
        help="path to directory with alternative Newick files, each with a single ultrametric tree",
    )

    parser.add_argument(
        "--num-to-select", default=2, type=int, help="the number of taxa to select"
    )
    parser.add_argument(
        "--use-patristic-distance-matrices",
        default=False,
        action="store_true",
        help="Developer option for switching between algorithms. "
        "Note that all input trees must be ultrametric, but some algorithms "
        "exploit this more effectively. Using these flags turns off these algorithms",
    )
    parser.add_argument(
        "--ultrametricity-tol",
        default=5e-5,
        type=float,
        help="precision for checking that the trees are ultrametric. Set this higher if you think that rounding error is causing your trees to be rejected as non-ultrametric",
    )
    parser.add_argument(
        "--max-solver-seconds",
        default=6000,
        type=float,
        help="Number of seconds that the solver is allowed to work on a single problem of maximizing the groupings over different trees (tree-dir mode only)",
    )
    parser.add_argument(
        "--scratch-dir",
        default=None,
        required=False,
        help="Directory from a previous run that was aborted.",
    )
    args = parser.parse_args(sys.argv[1:])
    if args.name_mapping_file is None:
        if args.country_file is not None:
            sys.exit("--country-file and --name-mapping-file must be used together.\n")
    else:
        if args.country_file is None:
            sys.exit("--country-file and --name-mapping-file must be used together.\n")
    if (args.tree_file is None) and (args.tree_dir is None):
        sys.exit("Either --tree-file or --tree-dir must be supplied.\n")
    if (args.tree_file is not None) and (args.tree_dir is not None):
        sys.exit("Only 1 of --tree-file or --tree-dir can be supplied.\n")
    if args.ultrametricity_tol < 0.0:
        sys.exit("--ultrametricity-tol cannot be negative")
    rs = RunSettings(
        country_name_fp=args.country_file,
        centroid_fp=args.centroid_file,
        tree_fp=args.tree_file,
        name_mapping_fp=args.name_mapping_file,
        num_to_select=args.num_to_select,
        use_ultrametricity=not args.use_patristic_distance_matrices,
        clade_defs_fp=args.clade_defs_file,
        name_updating_fp=args.name_updating_file,
        cut_branches_fp=args.cut_branches_file,
        tree_dir=args.tree_dir,
        ultrametric_tol=args.ultrametricity_tol,
        scratch_dir=args.scratch_dir,
        max_solver_seconds=args.max_solver_seconds,
    )
    return run(rs)


if __name__ == "__main__":
    rc = main()
    sys.exit(rc)
