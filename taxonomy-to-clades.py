#!/usr/bin/env pyth
import sys
import csv
from enum import Enum


class Ranks(Enum):
    SUBCLASS = 0
    INFRACLASS = 1
    MAGNORDER = 2
    SUPERORDER = 3
    ORDER = 4
    SUBORDER = 5
    INFRAORDER = 6
    PARVORDER = 7
    SUPERFAMILY = 8
    FAMILY = 9
    SUBFAMILY = 10
    TRIBE = 11
    GENUS = 12
    SUBGENUS = 13


def read_taxonomy_stream(inp):
    headers = [
        "sciName",
        "id",
        "phylosort",
        "mainCommonName",
        "otherCommonNames",
        "subclass",
        "infraclass",
        "magnorder",
        "superorder",
        "order",
        "suborder",
        "infraorder",
        "parvorder",
        "superfamily",
        "family",
        "subfamily",
        "tribe",
        "genus",
        "subgenus",
        "specificEpithet",
        "authoritySpeciesAuthor",
        "authoritySpeciesYear",
        "authorityParentheses",
        "originalNameCombination",
        "authoritySpeciesCitation",
        "authoritySpeciesLink",
        "holotypeVoucher",
        "holotypeVoucherURIs",
        "typeLocality",
        "typeLocalityLatitude",
        "typeLocalityLongitude",
        "nominalNames",
        "taxonomyNotes",
        "taxonomyNotesCitation",
        "distributionNotes",
        "distributionNotesCitation",
        "subregionDistribution",
        "countryDistribution",
        "continentDistribution",
        "biogeographicRealm",
        "iucnStatus",
        "extinct",
        "domestic",
        "flagged",
        "CMW_sciName",
        "diffSinceCMW",
        "MSW3_matchtype",
        "MSW3_sciName",
        "diffSinceMSW3",
    ]
    reader = csv.reader(inp, delimiter=",")
    clades_by_rank = [dict() for i in range(14)]
    for n, row in enumerate(reader):
        if n == 0:
            while row[-1] == "":
                del row[-1]
            if row != headers:
                raise RuntimeError(
                    "Expecting MDD headers:\n{}\n".format("\n".join(headers))
                )
            continue
        _process_row_into_cbr(row, clades_by_rank)
    print(Ranks.SUBCLASS.name, clades_by_rank[Ranks.SUBCLASS.value])
    print(Ranks.SUBCLASS.name, clades_by_rank[Ranks.SUBCLASS.value])
    print(Ranks.MAGNORDER.name, clades_by_rank[Ranks.MAGNORDER.value])
    print(Ranks.SUPERORDER.name, clades_by_rank[Ranks.SUPERORDER.value])
    print(Ranks.ORDER.name, clades_by_rank[Ranks.ORDER.value])
    print(Ranks.SUBORDER.name, clades_by_rank[Ranks.SUBORDER.value])
    print(Ranks.INFRAORDER.name, clades_by_rank[Ranks.INFRAORDER.value])
    print(Ranks.PARVORDER.name, clades_by_rank[Ranks.PARVORDER.value])
    print(Ranks.SUPERFAMILY.name, clades_by_rank[Ranks.SUPERFAMILY.value])
    print(Ranks.FAMILY.name, clades_by_rank[Ranks.FAMILY.value])
    print(Ranks.SUBFAMILY.name, clades_by_rank[Ranks.SUBFAMILY.value])
    print(Ranks.TRIBE.name, clades_by_rank[Ranks.TRIBE.value])
    print(Ranks.GENUS.name, clades_by_rank[Ranks.GENUS.value])
    print(Ranks.SUBGENUS.name, clades_by_rank[Ranks.SUBGENUS.value])
    all_clades = {}
    for rank in Ranks:
        ncbr = clades_by_rank[rank.value]
        rk = set(ncbr.keys())
        acs = set(all_clades.keys())
        iset = rk.intersection(acs)
        if iset:
            sys.stderr.write(
                f"Repeated clade names found in {rank.name} and previous ranks:\n"
            )
            for i in iset:
                sys.stderr.write(f"  {i}\n")
            assert False
        all_clades.update(ncbr)


_unset = frozenset(["NA", "INCERTAE SEDIS"])


def is_null(value):
    return value in _unset


def _process_row_into_cbr(row, clades_by_rank):
    sci_name = row[0]
    subclass = row[5 + Ranks.SUBCLASS.value]
    infraclass = row[5 + Ranks.INFRACLASS.value]
    magnorder = row[5 + Ranks.MAGNORDER.value]
    superorder = row[5 + Ranks.SUPERORDER.value]
    order = row[5 + Ranks.ORDER.value]
    suborder = row[5 + Ranks.SUBORDER.value]
    infraorder = row[5 + Ranks.INFRAORDER.value]
    parvorder = row[5 + Ranks.PARVORDER.value]
    superfamily = row[5 + Ranks.SUPERFAMILY.value]
    family = row[5 + Ranks.FAMILY.value]
    subfamily = row[5 + Ranks.SUBFAMILY.value]
    tribe = row[5 + Ranks.TRIBE.value]
    genus = row[5 + Ranks.GENUS.value]
    subgenus = row[5 + Ranks.SUBGENUS.value]
    assert not is_null(subclass.upper())
    assert not is_null(infraclass.upper())
    assert not is_null(order.upper())
    assert not is_null(family.upper())
    assert not is_null(genus.upper())
    subclass_d = clades_by_rank[Ranks.SUBCLASS.value]
    subclass_d.setdefault(subclass, set()).add(infraclass)
    infraclass_d = clades_by_rank[Ranks.INFRACLASS.value]
    if is_null(magnorder):
        key4infra = order if is_null(superorder) else superorder
    else:
        key4infra = magnorder
    infraclass_d.setdefault(infraclass, set()).add(key4infra)
    if not is_null(magnorder):
        magnorder_d = clades_by_rank[Ranks.MAGNORDER.value]
        key4magna = order if is_null(superorder) else superorder
        magnorder_d.setdefault(magnorder, set()).add(key4magna)
    if not is_null(superorder):
        superorder_d = clades_by_rank[Ranks.SUPERORDER.value]
        superorder_d.setdefault(superorder, set()).add(order)
    order_d = clades_by_rank[Ranks.ORDER.value]
    if not is_null(suborder):
        key4order = suborder
    else:
        if not is_null(infraorder):
            key4order = infraorder
        else:
            if not is_null(parvorder):
                key4order = parvorder
            else:
                key4order = family if is_null(superfamily) else superfamily
    order_d.setdefault(order, set()).add(key4order)
    if not is_null(suborder):
        suborder_d = clades_by_rank[Ranks.SUBORDER.value]
        if not is_null(infraorder):
            key4suborder = infraorder
        else:
            if not is_null(parvorder):
                key4suborder = parvorder
            else:
                key4suborder = family if is_null(superfamily) else superfamily
        suborder_d.setdefault(suborder, set()).add(key4suborder)
    if not is_null(infraorder):
        infraorder_d = clades_by_rank[Ranks.INFRAORDER.value]
        if not is_null(parvorder):
            key4infraorder = parvorder
        else:
            key4infraorder = family if is_null(superfamily) else superfamily
        infraorder_d.setdefault(infraorder, set()).add(key4infraorder)
    if not is_null(parvorder):
        parorder_d = clades_by_rank[Ranks.PARVORDER.value]
        key4parvorder = family if is_null(superfamily) else superfamily
        parorder_d.setdefault(parvorder, set()).add(key4parvorder)
    if not is_null(superfamily):
        superfamily_d = clades_by_rank[Ranks.SUPERFAMILY.value]
        superfamily_d.setdefault(superfamily, set()).add(family)
    family_d = clades_by_rank[Ranks.FAMILY.value]
    if not is_null(subfamily):
        key4family = subfamily
    else:
        key4family = genus if is_null(tribe) else tribe
    family_d.setdefault(family, set()).add(key4family)
    if not is_null(subfamily):
        subfamily_d = clades_by_rank[Ranks.SUBFAMILY.value]
        key4subfamily = genus if is_null(tribe) else tribe
        subfamily_d.setdefault(subfamily, set()).add(key4subfamily)
    if not is_null(tribe):
        tribe_d = clades_by_rank[Ranks.TRIBE.value]
        tribe_d.setdefault(tribe, set()).add(genus)
    genus_d = clades_by_rank[Ranks.GENUS.value]
    key4genus = sci_name if is_null(subgenus) else f"{genus } subgenus {subgenus}"
    genus_d.setdefault(genus, set()).add(key4genus)
    if not is_null(subgenus):
        subgenus_d = clades_by_rank[Ranks.SUBGENUS.value]
        subgenus_d.setdefault(f"{genus } subgenus {subgenus}", set()).add(sci_name)


def main(taxonomy_fp=None):
    if taxonomy_fp is None:
        read_taxonomy_stream(sys.stdin)
    else:
        with open(taxonomy_fp, "r") as inp:
            read_taxonomy_stream(inp)


if __name__ == "__main__":
    main(None if len(sys.argv) == 1 else sys.argv[1])
