#! /usr/bin/env python3
import dendropy
import csv
from .logs import info
from .taxonomy import parse_clade_defs
from .tree_cleaning import prune_taxa_without_sp_data


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
                info(f'Skipping taxon "{sp_name}" due to NA in centroid.')
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


def parse_name_updating(name_updating_fp):
    mapping = {}
    rev_map_set = set()
    if not name_updating_fp:
        return {}
    with open(name_updating_fp, "r", newline="") as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for n, row in enumerate(reader):
            if n == 0:
                if row[0] != "CMW_sciName" or row[1] != "sciName":
                    raise RuntimeError(
                        f"Expecting {name_updating_fp} to have CMW_sciName and sciName as the headers"
                        "as produced by taxonomy-to-clades.py"
                    )
                continue
            tree_label, tax_label = row[0], row[1]
            if tree_label in mapping:
                raise RuntimeError(
                    f"{tree_label} repeated in first column of {name_updating_fp}"
                )
            if tax_label in rev_map_set:
                raise RuntimeError(
                    f"{tax_label} repeated in second column of {name_updating_fp}"
                )
            mapping[tree_label] = tax_label
            rev_map_set.add(tax_label)
    return mapping


def parse_geo_and_tree(
    country_name_fp,
    centroid_fp,
    name_mapping_fp,
    tree_fp,
    clade_defs_fp,
    name_updating_fp=None,
):
    new_names_for_leaves = parse_name_updating(name_updating_fp)
    clades = parse_clade_defs(clade_defs_fp)
    tree = dendropy.Tree.get(path=tree_fp, schema="nexus")
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

    prune_taxa_without_sp_data(
        tree,
        frozenset(sp_by_name.keys()),
        upham_to_iucn=upham_to_iucn,
        name_mapping_fp=name_mapping_fp,
        centroid_fp=centroid_fp,
        clades=clades,
        new_names_for_leaves=new_names_for_leaves,
    )
    return tree, sp_by_name
