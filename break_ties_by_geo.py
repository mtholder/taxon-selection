#!/usr/bin/env python
import sys

from geotaxsel import parse_geo, min_dist_between_sp


def error(msg):
    sys.stderr.write(f"break_ties_by_geo.py: ERROR: {msg}\n")


def status(msg):
    sys.stderr.write(f"{msg}\n")


def read_chosen_taxa(chosen_tax_fp):
    names = []
    with open(chosen_tax_fp, "r") as inp:
        names = [i.strip() for i in inp if i.strip()]
    return names


def parse_cut_branches_file(cut_branches_fp):
    sets_by_members = {}
    num_cuts = 0
    with open(cut_branches_fp, "r") as inp:
        header_read = False
        for line in inp:
            if header_read is False:
                assert line == "num-leaves\tnames\n"
                header_read = True
                continue
            if not line.strip():
                continue
            bl = line.strip().split("\t")
            assert len(bl) == 2
            num_members = int(bl[0])
            choice_set = frozenset([i.strip() for i in bl[1].split(",")])
            assert len(choice_set) == num_members
            for member in choice_set:
                if member in sets_by_members:
                    raise RuntimeError(
                        f"multiple sets in cut branches have taxon {member}"
                    )
                sets_by_members[member] = choice_set
            num_cuts += 1
    return sets_by_members, num_cuts


def main(cut_branches_fp, chosen_tax_fp, centroid_fp):
    try:
        blob = parse_geo(
            country_name_fp=None,
            centroid_fp=centroid_fp,
            name_mapping_fp=None,
            clade_defs_fp=None,
            name_updating_fp=None,
        )
    except:
        error(f"Problem reading centroid_fp {centroid_fp}")
        raise
    sp_by_name = blob[0]
    chosen = read_chosen_taxa(chosen_tax_fp)
    status(f"{len(chosen)} chosen taxa.")
    sp_to_set, num_cuts = parse_cut_branches_file(cut_branches_fp)
    status(f"{num_cuts} cuts read.")
    assert num_cuts == len(chosen)
    to_print = []
    max_num_choices = 0
    for first_choice in chosen:
        assert first_choice in sp_to_set
        choice_set = sp_to_set[first_choice]
        dist_ch_list = [(-1, first_choice)]
        fsp = sp_by_name[first_choice]
        for el in choice_set:
            if el != first_choice:
                el_sp = sp_by_name[el]
                md = min_dist_between_sp(fsp, el_sp)
                dist_ch_list.append((md, el))
        dist_ch_list.sort()
        ranked = [i[1] for i in dist_ch_list]
        assert ranked[0] == first_choice
        assert len(ranked) == len(choice_set)
        if len(ranked) > max_num_choices:
            max_num_choices = len(ranked)
        to_print.append("\t".join(ranked))
    header = [f"choice-{i}" for i in range(1, 1 + max_num_choices)]
    print("\t".join(header))
    to_print.sort()
    for line in to_print:
        print(line)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("Expecting 3 arguments: cut_branches_fp chosen_tax_fp centroid_fp")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
