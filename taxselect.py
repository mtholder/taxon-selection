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
    choose_most_common,
)
import dendropy


def record_clade_sel(sel, rep_selections):
    for anc in sel:
        leaves_below = list(anc.leaf_nodes())
        labels_below = frozenset([i.taxon.label for i in leaves_below])
        pn = rep_selections.get(labels_below, 0)
        rep_selections[labels_below] = 1 + pn


def run_tree_dir(
    country_name_fp,
    centroid_fp,
    name_mapping_fp,
    num_to_select,
    use_ultrametricity=True,
    clade_defs_fp=None,
    name_updating_fp=None,
    cut_branches_fp=None,
    tree_dir=None,
    ultrametric_tol=5e-5,
):
    geo_ret = parse_geo(
        country_name_fp=country_name_fp,
        centroid_fp=centroid_fp,
        name_mapping_fp=name_mapping_fp,
        clade_defs_fp=clade_defs_fp,
        name_updating_fp=name_updating_fp,
    )
    sp_by_name, clades, upham_to_iucn, new_names_for_leaves = geo_ret
    file_names = os.listdir(tree_dir)
    file_names.sort()
    sp_pat = re.compile("^([A-Z][a-z]+ +[-a-z0-9]+)$")
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
        if tree_n > 9:
            break
    final = choose_most_common(rep_selections, num_to_select, num_trees=1 + tree_n)
    output_chosen_anc(
        tree=None,
        cut_branches_fp=cut_branches_fp,
        chosen_ancs=final,
        label_sets_to_freq=rep_selections,
    )
    return


def run(
    country_name_fp,
    centroid_fp,
    name_mapping_fp,
    num_to_select,
    tree_fp=None,
    use_ultrametricity=True,
    clade_defs_fp=None,
    name_updating_fp=None,
    cut_branches_fp=None,
    tree_dir=None,
    ultrametric_tol=5e-5,
):
    if tree_dir is not None:
        return run_tree_dir(
            country_name_fp=country_name_fp,
            centroid_fp=centroid_fp,
            name_mapping_fp=name_mapping_fp,
            num_to_select=num_to_select,
            use_ultrametricity=use_ultrametricity,
            clade_defs_fp=clade_defs_fp,
            name_updating_fp=name_updating_fp,
            cut_branches_fp=cut_branches_fp,
            tree_dir=tree_dir,
            ultrametric_tol=ultrametric_tol,
        )
    sp_pat = re.compile("^([A-Z][a-z]+ +[-a-z0-9]+) [A-Z][A-Za-z]+ [A-Z]+$")
    assert tree_fp is not None
    tree, sp_by_name = parse_geo_and_tree(
        country_name_fp,
        centroid_fp,
        name_mapping_fp,
        tree_fp,
        clade_defs_fp=clade_defs_fp,
        name_updating_fp=name_updating_fp,
        sp_pat_in_tree=sp_pat,
    )
    if use_ultrametricity:
        sel = ultrametric_greedy_mmd(
            tree, num_to_select, sp_by_name, ultrametric_tol=ultrametric_tol
        )
    else:
        sel = greedy_mmd(tree, num_to_select, sp_by_name)
    output_chosen_anc(tree, cut_branches_fp, sel)
    sys.exit("early exit\n")
    print("Selected:\n  {}\n".format("\n  ".join(sel)))
    return


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
    run(
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
    )


if __name__ == "__main__":
    main()
