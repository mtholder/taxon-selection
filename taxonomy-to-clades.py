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
    assert subclass != "NA"
    assert infraclass != "NA"
    subclass_d = clades_by_rank[Ranks.SUBCLASS.value]
    subclass_d.setdefault(subclass, set()).add(infraclass)
    infraclass_d = clades_by_rank[Ranks.INFRACLASS.value]
    assert order != "NA"
    if magnorder == "NA":
        if superorder == "NA":
            key4infra = order
        else:
            key4infra = superorder
    else:
        key4infra = magnorder
    infraclass_d.setdefault(infraclass, set()).add(key4infra)
    if magnorder != "NA":
        magnorder_d = clades_by_rank[Ranks.MAGNORDER.value]
        if superorder == "NA":
            key4magna = order
        else:
            key4magna = superorder
        magnorder_d.setdefault(magnorder, set()).add(key4magna)
    if superorder != "NA":
        superorder_d = clades_by_rank[Ranks.SUPERORDER.value]
        superorder_d.setdefault(superorder, set()).add(order)


def main(taxonomy_fp=None):
    if taxonomy_fp is None:
        read_taxonomy_stream(sys.stdin)
    else:
        with open(taxonomy_fp, "r") as inp:
            read_taxonomy_stream(inp)


if __name__ == "__main__":
    main(None if len(sys.argv) == 1 else sys.argv[1])
