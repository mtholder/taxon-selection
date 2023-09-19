#!/usr/bin/env python
import dendropy
import csv
import sys
import re
from geopy.distance import geodesic
from dendropy.calculate.phylogeneticdistance import PhylogeneticDistanceMatrix


def calc_dist(loc_1, loc_2):
    if loc_1 is loc_2:
        return 0.0
    return geodesic(loc_1.coords, loc_2.coords).km


def read_country_names(country_name_fp):
    countries = []
    with open(country_name_fp, "r", newline="", encoding="latin-1") as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        for n, row in enumerate(reader):
            if n == 0:
                assert row[0] == "COUNTRY"
                continue
            countries.append(row[0].strip())
    return countries


class Loc(object):
    def __init__(self, latitude, longitude):
        self.str_lat = latitude
        self.str_long = longitude
        self.latitude = float(latitude)
        self.longitude = float(longitude)
        self.coords = (self.latitude, self.longitude)

    def __str__(self):
        return "Loc({}, {})".format(self.str_lat, self.str_long)

    def __hash__(self):
        return hash((self.str_lat, self.str_long))


class NamedLoc(object):
    def __init__(self, name, c_id, loc):
        self.name = name
        self.id = c_id
        self.loc = loc
        self._h = hash((self.name, self.id, self.loc))
        self.idx = None
        self.coords = loc.coords
        # _all_loc.append(self)

    def __str__(self):
        return "Country({}, {}, {})".format(self.name, self.id, self.loc)

    def __hash__(self):
        return self._h


class Species(object):
    def __init__(self, name, sp_id):
        self.name = name
        self.id = sp_id
        self.locations = set()

    def add_loc(self, loc):
        self.locations.add(loc)


def read_centroids(centroid_fp, countries):
    if countries is None:
        return read_centroids_sans_countries(centroid_fp)
    return read_centroids_with_countries(centroid_fp, countries)


def read_centroids_sans_countries(centroid_fp):
    sp_by_name = {}
    with open(centroid_fp, "r", newline="", encoding="latin-1") as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        for n, row in enumerate(reader):
            if n == 0:
                assert row == ["Upham_name", "x", "y"]
                continue
            sp_name, longitude, latitude = row
            if longitude == "NA" or latitude == "NA":
                sys.stderr.write(f'Skipping taxon "{sp_name}" due to NA in centroid.\n')
                continue
            loc = Loc(latitude, longitude)
            this_loc = NamedLoc(name=None, c_id=None, loc=loc)
            sp = sp_by_name.get(sp_name)
            assert sp is None
            sp = Species(sp_name, sp_id=None)
            sp_by_name[sp_name] = sp
            sp.add_loc(this_loc)
    return sp_by_name


def read_centroids_with_countries(centroid_fp, countries):
    sp_by_name = {}
    country_by_name = {}
    with open(centroid_fp, "r", newline="", encoding="latin-1") as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        for n, row in enumerate(reader):
            if n == 0:
                assert row == [
                    "Species_no",
                    "binomial",
                    "Country_ID",
                    "Country",
                    "Longitude",
                    "Latitude",
                ]
                continue
            sp_n, sp_name, countr_id, countr_name, longitude, latitude = row
            assert countr_name in countries
            pc_list = country_by_name.setdefault(countr_name, [])
            this_loc = None
            if pc_list is not None:
                for pc in pc_list:
                    if (
                        pc.name == countr_name
                        and pc.id == countr_id
                        and pc.loc.str_lat == latitude
                        and pc.loc.str_long == longitude
                    ):
                        this_loc = pc
                        break
            if this_loc is None:
                loc = Loc(latitude, longitude)
                this_loc = NamedLoc(countr_name, countr_id, loc)
                pc_list.append(this_loc)
            sp = sp_by_name.get(sp_name)
            if sp is None:
                sp = Species(sp_name, sp_n)
                sp_by_name[sp_name] = sp
            sp.add_loc(this_loc)
    return sp_by_name


def read_upham_to_iucn(name_mapping_fp):
    up_to_iucn = {}
    with open(name_mapping_fp, "r", newline="", encoding="utf-8") as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for n, row in enumerate(reader):
            if n == 0:
                assert [i.lower() for i in row] == [
                    "upham_name",
                    "iucn_name",
                    "total_avail",
                ]
                continue
            upham, iucn, tot_avail = row
            upham = " ".join(upham.split("_"))
            if upham == "NA":
                continue
            iucn = " ".join(iucn.split("_"))
            if upham in up_to_iucn:
                raise RuntimeError("repeated name '{}'".format(row[0]))
            up_to_iucn[upham] = iucn
    return up_to_iucn


VERBOSE = True


def debug(m):
    if VERBOSE:
        sys.stderr.write(f"greedy_mmd debug: {m}\n")


def sel_most_geo_div_taxon(label_ind_pairs, loc_list, sp_by_name):
    # sys.stderr.write('passed in loc_list = {}\n'.format(loc_list))
    md = None
    md_loc = None
    md_sp = None
    for label, ind in label_ind_pairs:
        sp = sp_by_name[label]
        for l1 in sp.locations:
            sum_sq_dist = 0.0
            for l2 in loc_list:
                try:
                    sum_sq_dist += calc_dist(l1, l2)
                except:
                    # sys.stderr.write('sp.__dict__ = {}\n'.format(sp.__dict__))
                    # sys.stderr.write('loc_list = {}\n'.format(loc_list))
                    # sys.stderr.write('l1 = {}\n'.format(l1))
                    # sys.stderr.write('l2 = {}\n'.format(l2))
                    raise
            if md is None or sum_sq_dist > md:
                md = sum_sq_dist
                md_loc = l1
                md_sp_ind = ind
    assert md_sp_ind is not None
    return md_sp_ind, md_loc


def most_divergent_locs(tax_1, tax_2, sp_by_name):
    sp1, sp2 = sp_by_name[tax_1], sp_by_name[tax_2]
    md = None
    md_pair = None
    for loc1 in sp1.locations:
        for loc2 in sp2.locations:
            d = calc_dist(loc1, loc2)
            if md is None or d > md:
                md = d
                md_pair = (loc1, loc2)
    return md_pair


def greedy_mmd(tree, num_taxa, sp_by_name):
    # not efficiently updating min dist yet
    taxa_list = [i.taxon for i in tree.leaf_nodes()]
    taxa_label_list = [i.label for i in taxa_list]
    if num_taxa == len(taxa_label_list):
        return taxa_label_list
    if num_taxa > len(taxa_label_list):
        raise ValueError("num_taxa exceeds the number of taxa in the tree")
    debug(f"Calculating patristic distance matrix")
    pdmc = PhylogeneticDistanceMatrix.from_tree(tree)
    pdm = []
    max_dist = None
    max_dist_inds = None
    num_all_tax = len(taxa_list)
    debug(f"Copying patristic distance matrix for {num_all_tax}")
    for row_ind, row_tax in enumerate(taxa_list):
        row = []
        if row_ind % 100 == 0:
            debug(f"  row {row_ind}")
        for col_ind, col_tax in enumerate(taxa_list):
            nd = pdmc.patristic_distance(row_tax, col_tax)
            if max_dist is None or nd > max_dist:
                max_dist = nd
                max_dist_inds = (row_ind, col_ind)
            row.append(nd)
        pdm.append(row)
    sel_tax_labels = [
        taxa_label_list[max_dist_inds[0]],
        taxa_label_list[max_dist_inds[1]],
    ]
    sel_inds = set(max_dist_inds)
    debug(f'Most divergent 2 taxa are "{sel_tax_labels}" with dist= {max_dist}')
    tax_1, tax_2 = taxa_list[max_dist_inds[0]], taxa_list[max_dist_inds[1]]
    loc1, loc2 = most_divergent_locs(tax_1.label, tax_2.label, sp_by_name)
    TOL = 1.0e-5
    curr_locs = [loc1, loc2]
    while len(sel_inds) < num_taxa:
        debug(f"Finding taxon {1 + len(sel_inds)}...")
        max_min_dist, mmd_ind = -1, None
        mmd_ind_set = set()
        for row_ind, row_tax in enumerate(taxa_list):
            if row_ind in sel_inds:
                continue
            row = pdm[row_ind]
            inc_dist = [row[i] for i in sel_inds]
            curr_min_dist = min(inc_dist)
            if curr_min_dist > max_min_dist:
                max_min_dist = curr_min_dist
                mmd_ind = row_ind
                mmd_ind_set = set()
                mmd_ind_set.add(row_ind)
            elif abs(curr_min_dist - max_min_dist) < TOL:
                mmd_ind_set.add(row_ind)
        tied_tax = [(taxa_label_list[i], i) for i in mmd_ind_set]
        mmd_ind, sel_loc = sel_most_geo_div_taxon(tied_tax, curr_locs, sp_by_name)
        ntl = taxa_label_list[mmd_ind]
        debug(
            f'    taxon {1 + len(sel_inds)} = "{ntl}" with MD = {curr_min_dist} set = {mmd_ind_set}'
        )
        sel_tax_labels.append(ntl)
        sel_inds.add(mmd_ind)
        curr_locs.append(sel_loc)
    debug("Taxa selected, cleaning up...")
    return sel_tax_labels


def main(country_name_fp, centroid_fp, name_mapping_fp, tree_fp, num_to_select):
    if country_name_fp is not None:
        countries = read_country_names(country_name_fp)
        countries = frozenset(countries)
        assert name_mapping_fp is not None
        upham_to_iucn = read_upham_to_iucn(name_mapping_fp)
    else:
        assert name_mapping_fp is None
        countries = None
        upham_to_iucn = None
    sp_by_name = read_centroids(centroid_fp, countries)
    tree = dendropy.Tree.get(path=tree_fp, schema="nexus")
    taxa_list = [i for i in tree.taxon_namespace]
    to_prune = []
    sp_pat = re.compile("^([A-Z][a-z]+ +[-a-z0-9]+) [A-Z][A-Za-z]+ [A-Z]+$")
    for n, i in enumerate(taxa_list):
        m = sp_pat.match(i.label)
        if not m:
            to_prune.append(i)
            sys.stderr.write(
                'Will prune "{}". did not match species name pattern\n'.format(
                    i.label, name_mapping_fp
                )
            )
            continue
        sp_name = m.group(1)
        if upham_to_iucn is not None:
            final_name = upham_to_iucn.get(sp_name)
        else:
            final_name = sp_name
        if final_name is None:
            to_prune.append(i)
            sys.stderr.write(
                'Will prune "{}" not found in {}\n'.format(sp_name, name_mapping_fp)
            )
        elif final_name not in sp_by_name:
            to_prune.append(i)
            sys.stderr.write(
                'Will prune "{}" not found in {}\n'.format(sp_name, centroid_fp)
            )
        else:
            i.label = final_name
    tree.prune_taxa(to_prune)
    sel = greedy_mmd(tree, num_to_select, sp_by_name)
    print("Selected:\n  {}\n".format("\n  ".join(sel)))
    return


if __name__ == "__main__":
    if len(sys.argv) == 6:
        main(
            country_name_fp=sys.argv[1],
            centroid_fp=sys.argv[2],
            name_mapping_fp=sys.argv[3],
            tree_fp=sys.argv[4],
            num_to_select=int(sys.argv[5]),
        )
    else:
        assert len(sys.argv) == 4
        main(
            country_name_fp=None,
            centroid_fp=sys.argv[1],
            tree_fp=sys.argv[2],
            name_mapping_fp=None,
            num_to_select=int(sys.argv[3]),
        )
