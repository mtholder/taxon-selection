#!/usr/bin/env python
import argparse
import sys
from geotaxsel import greedy_mmd, parse_geo_and_tree, ultrametric_greedy_mmd


def run(
    country_name_fp,
    centroid_fp,
    name_mapping_fp,
    tree_fp,
    num_to_select,
    use_ultrametricity=True,
    clade_defs_fp=None,
    name_updating_fp=None,
    cut_branches_fp=None,
):
    tree, sp_by_name = parse_geo_and_tree(
        country_name_fp,
        centroid_fp,
        name_mapping_fp,
        tree_fp,
        clade_defs_fp=clade_defs_fp,
        name_updating_fp=name_updating_fp,
    )
    if use_ultrametricity:
        sel = ultrametric_greedy_mmd(tree, num_to_select, sp_by_name, cut_branches_fp)
    else:
        sel = greedy_mmd(tree, num_to_select, sp_by_name)
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
        required=True,
        help="path to NEXUS file with a single ultrametric tree",
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
    args = parser.parse_args(sys.argv[1:])
    if args.name_mapping_file is None:
        if args.country_file is not None:
            sys.exit("--country-file and --name-mapping-file must be used together.\n")
    else:
        if args.country_file is None:
            sys.exit("--country-file and --name-mapping-file must be used together.\n")
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
    )


if __name__ == "__main__":
    main()
