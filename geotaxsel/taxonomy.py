#! /usr/bin/env python3
import csv
import sys
from enum import Enum
import re


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


_rank_str2rank_dict = {_i.name: _i for _i in Ranks}


def rank_str_to_rank(s):
    if s is None:
        return None
    try:
        return _rank_str2rank_dict[s]
    except KeyError:
        raise ValueError(f"Rank {repr(s)} not recognized.")


class CladeDef(object):
    def __init__(self, mrca_of=None, incertae_sedis=None, rank=None):
        self.must = set()
        self.might = set()
        self.rank = rank
        self.name = None
        if mrca_of:
            self.must.update(mrca_of)
        if incertae_sedis:
            self.might.update(incertae_sedis)

    def add(self, element):
        self.must.add(element)

    def add_possible(self, element):
        self.might.add(element)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        rs = ""
        if self.rank is not None:
            rs = f", rank=Ranks.{self.rank.name}"
        if self.might:
            return f"CladeDef(mrca_of={repr(self.must)}, incertae_sedis={repr(self.might)}{rs})"
        return f"CladeDef(mrca_of={repr(self.must)}{rs})"


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
    # noinspection PyUnusedLocal
    clades_by_rank = [dict() for i in range(14)]
    incertae_sedis = []
    sci_name_idx = headers.index("sciName")
    cmw_sci_name_idx = headers.index("CMW_sciName")
    cmw2current = {}
    for n, row in enumerate(reader):
        if n == 0:
            while row[-1] == "":
                del row[-1]
            if row != headers:
                raise RuntimeError(
                    "Expecting MDD headers:\n{}\n".format("\n".join(headers))
                )
            continue
        raw_sci_name = row[sci_name_idx].strip()
        raw_cmw_name = row[cmw_sci_name_idx].strip()
        clean_sci = " ".join(raw_sci_name.split("_"))
        clean_cmw = " ".join(raw_cmw_name.split("_"))
        if clean_cmw and (clean_sci != clean_cmw):
            cmw2current[clean_cmw] = clean_sci
        _process_row_into_cbr(row, clades_by_rank, incertae_sedis)
    # print(Ranks.SUBCLASS.name, clades_by_rank[Ranks.SUBCLASS.value])
    # print(Ranks.SUBCLASS.name, clades_by_rank[Ranks.SUBCLASS.value])
    # print(Ranks.MAGNORDER.name, clades_by_rank[Ranks.MAGNORDER.value])
    # print(Ranks.SUPERORDER.name, clades_by_rank[Ranks.SUPERORDER.value])
    # print(Ranks.ORDER.name, clades_by_rank[Ranks.ORDER.value])
    # print(Ranks.SUBORDER.name, clades_by_rank[Ranks.SUBORDER.value])
    # print(Ranks.INFRAORDER.name, clades_by_rank[Ranks.INFRAORDER.value])
    # print(Ranks.PARVORDER.name, clades_by_rank[Ranks.PARVORDER.value])
    # print(Ranks.SUPERFAMILY.name, clades_by_rank[Ranks.SUPERFAMILY.value])
    # print(Ranks.FAMILY.name, clades_by_rank[Ranks.FAMILY.value])
    # print(Ranks.SUBFAMILY.name, clades_by_rank[Ranks.SUBFAMILY.value])
    # print(Ranks.TRIBE.name, clades_by_rank[Ranks.TRIBE.value])
    # print(Ranks.GENUS.name, clades_by_rank[Ranks.GENUS.value])
    # print(Ranks.SUBGENUS.name, clades_by_rank[Ranks.SUBGENUS.value])
    all_clades = {}
    for rank in Ranks:
        ncbr = clades_by_rank[rank.value]
        for v in ncbr.values():
            v.rank = rank
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
    for datum in incertae_sedis:
        _process_incertae_sedis(datum, all_clades)
    retlist = []
    for rank in Ranks:
        cbr = clades_by_rank[rank.value]
        skl = list(cbr.keys())
        skl.sort()
        for k in skl:
            v = cbr[k]
            retlist.append((k, v))
    return retlist, cmw2current


_unset = frozenset(["NA", "INCERTAE SEDIS"])


def is_null(value):
    return value.upper() in _unset


def is_inc_sedis(name):
    return name.upper() == "INCERTAE SEDIS"


def _process_incertae_sedis(data_row, clades_dict):
    (
        sci_name,
        subclass,
        infraclass,
        magnorder,
        superorder,
        order,
        suborder,
        infraorder,
        parvorder,
        superfamily,
        family,
        subfamily,
        tribe,
        genus,
        subgenus,
    ) = data_row
    # To make life simple we only deal with the cases we need to.
    assert not is_inc_sedis(subclass)
    assert not is_inc_sedis(infraclass)
    assert not is_inc_sedis(magnorder)
    assert not is_inc_sedis(superorder)
    assert not is_inc_sedis(order)
    assert not is_inc_sedis(suborder)
    assert not is_inc_sedis(infraorder)
    assert not is_inc_sedis(parvorder)
    assert not is_inc_sedis(family)
    assert not is_inc_sedis(genus)
    if is_inc_sedis(superfamily):
        if not is_null(parvorder):
            key4superis = parvorder
        else:
            assert not is_null(infraorder)
            key4superis = infraorder
        assert not is_null(family)
        clade_def = clades_dict[key4superis]
        for subc in clade_def.must:
            sc = clades_dict[subc]
            sc.add_possible(family)
    if is_inc_sedis(subfamily):
        assert not is_null(family)
        clade_def = clades_dict[family]
        if not is_null(tribe):
            val4fam = tribe
        else:
            assert not is_null(genus)
            val4fam = genus
        for subc in clade_def.must:
            sc = clades_dict[subc]
            sc.add_possible(val4fam)
    if is_inc_sedis(subgenus):
        assert not is_null(genus)
        clade_def = clades_dict[genus]
        for subc in clade_def.must:
            try:
                sc = clades_dict[subc]
            except KeyError:
                pass
            else:
                sc.add_possible(sci_name)
    if is_inc_sedis(subfamily):
        assert not is_null(family)
        clade_def = clades_dict[family]
        if not is_null(tribe):
            val4fam = tribe
        else:
            assert not is_null(genus)
            val4fam = genus
        for subc in clade_def.must:
            sc = clades_dict[subc]
            sc.add_possible(val4fam)


def _process_row_into_cbr(row, clades_by_rank, incertae_sedis_rows):
    sci_name = " ".join(row[0].split("_"))
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
    data = [
        sci_name,
        subclass,
        infraclass,
        magnorder,
        superorder,
        order,
        suborder,
        infraorder,
        parvorder,
        superfamily,
        family,
        subfamily,
        tribe,
        genus,
        subgenus,
    ]
    if any([i.upper() == "INCERTAE SEDIS" for i in data]):
        incertae_sedis_rows.append(data)
    assert not is_null(subclass.upper())
    assert not is_null(infraclass.upper())
    assert not is_null(order.upper())
    assert not is_null(family.upper())
    assert not is_null(genus.upper())
    # Register Subclass clade
    subclass_d = clades_by_rank[Ranks.SUBCLASS.value]
    subclass_d.setdefault(subclass, CladeDef()).add(infraclass)

    # Register Infraclass clade
    infraclass_d = clades_by_rank[Ranks.INFRACLASS.value]
    if not is_null(magnorder):
        key4infra = magnorder
    else:
        key4infra = order if is_null(superorder) else superorder
    infraclass_d.setdefault(infraclass, CladeDef()).add(key4infra)
    # Register Magnaorder clade
    if not is_null(magnorder):
        magnorder_d = clades_by_rank[Ranks.MAGNORDER.value]
        key4magna = order if is_null(superorder) else superorder
        magnorder_d.setdefault(magnorder, CladeDef()).add(key4magna)
    # Register Superorder clade
    if not is_null(superorder):
        superorder_d = clades_by_rank[Ranks.SUPERORDER.value]
        superorder_d.setdefault(superorder, CladeDef()).add(order)
    # Register Order clade
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
    order_d.setdefault(order, CladeDef()).add(key4order)
    # Register Suborder clade
    if not is_null(suborder):
        suborder_d = clades_by_rank[Ranks.SUBORDER.value]
        if not is_null(infraorder):
            key4suborder = infraorder
        else:
            if not is_null(parvorder):
                key4suborder = parvorder
            else:
                key4suborder = family if is_null(superfamily) else superfamily
        suborder_d.setdefault(suborder, CladeDef()).add(key4suborder)
    # Register Infraorder clade
    if not is_null(infraorder):
        infraorder_d = clades_by_rank[Ranks.INFRAORDER.value]
        if not is_null(parvorder):
            key4infraorder = parvorder
        else:
            key4infraorder = family if is_null(superfamily) else superfamily
        infraorder_d.setdefault(infraorder, CladeDef()).add(key4infraorder)
    # Register Parvorder clade
    if not is_null(parvorder):
        parorder_d = clades_by_rank[Ranks.PARVORDER.value]
        key4parvorder = family if is_null(superfamily) else superfamily
        parorder_d.setdefault(parvorder, CladeDef()).add(key4parvorder)
    # Register Superfamily clade
    if not is_null(superfamily):
        superfamily_d = clades_by_rank[Ranks.SUPERFAMILY.value]
        superfamily_d.setdefault(superfamily, CladeDef()).add(family)
    # Register Family clade
    family_d = clades_by_rank[Ranks.FAMILY.value]
    if not is_null(subfamily):
        key4family = subfamily
    else:
        key4family = genus if is_null(tribe) else tribe
    family_d.setdefault(family, CladeDef()).add(key4family)
    # Register Subfamily clade
    if not is_null(subfamily):
        subfamily_d = clades_by_rank[Ranks.SUBFAMILY.value]
        key4subfamily = genus if is_null(tribe) else tribe
        subfamily_d.setdefault(subfamily, CladeDef()).add(key4subfamily)
    # Register Subfamily clade
    if not is_null(tribe):
        tribe_d = clades_by_rank[Ranks.TRIBE.value]
        tribe_d.setdefault(tribe, CladeDef()).add(genus)
    # Register Genus clade
    genus_d = clades_by_rank[Ranks.GENUS.value]
    key4genus = sci_name if is_null(subgenus) else f"{genus } subgenus {subgenus}"
    genus_d.setdefault(genus, CladeDef()).add(key4genus)
    # Register Subgenus clade
    if not is_null(subgenus):
        subgenus_d = clades_by_rank[Ranks.SUBGENUS.value]
        subgenus_d.setdefault(f"{genus } subgenus {subgenus}", CladeDef()).add(sci_name)


_pref_pat = r"^\s*CladeDef\(\s*mrca_of\s*=\s*\{([^}]+)\}\s*"
_inc_sed_pat = r",\s*incertae_sedis\s*=\s*\{([^}]+)\}\s*"
_ranks_pat = r",\s*rank\s*=\s*Ranks[.]([-_A-Za-z]+)\s*"
_close_parens_pat = r"\)\s*$"
_cdef_arg_pat = re.compile(f"{_pref_pat}{_close_parens_pat}")
_cdef_is_arg_pat = re.compile(f"{_pref_pat}{_inc_sed_pat}{_close_parens_pat}")
_cdef_r_arg_pat = re.compile(f"{_pref_pat}{_ranks_pat}{_close_parens_pat}")
_cdef_is_r_arg_pat = re.compile(
    f"{_pref_pat}{_inc_sed_pat}{_ranks_pat}{_close_parens_pat}"
)


def _parse_cdef_str(defn):
    m = _cdef_r_arg_pat.match(defn)
    if m:
        mrca = td_set_str_contents_to_list(m.group(1))
        rank = rank_str_to_rank(m.group(2))
        return CladeDef(mrca_of=mrca, rank=rank)
    m = _cdef_is_r_arg_pat.match(defn)
    if m:
        mrca = td_set_str_contents_to_list(m.group(1))
        inc_sed = td_set_str_contents_to_list(m.group(2))
        rank = rank_str_to_rank(m.group(3))
        return CladeDef(mrca_of=mrca, incertae_sedis=inc_sed, rank=rank)
    m = _cdef_arg_pat.match(defn)
    if m:
        mrca = td_set_str_contents_to_list(m.group(1))
        return CladeDef(mrca_of=mrca)
    m = _cdef_is_arg_pat.match(defn)
    if m:
        mrca = td_set_str_contents_to_list(m.group(1))
        inc_sed = td_set_str_contents_to_list(m.group(2))
        return CladeDef(mrca_of=mrca, incertae_sedis=inc_sed)
    raise RuntimeError(f"Could not parse: {defn}")


def td_set_str_contents_to_list(str_contents):
    inp_list = str_contents.split(",")
    out_list = []
    for unclean in inp_list:
        word = unclean.strip()
        ok = False
        if word[0] == "'":
            ok = word[-1] == "'"
        elif word[0] == '"':
            ok = word[-1] != '"'
        if not ok:
            raise RuntimeError(f"{word} not quoted easily.")
        unquoted = word[1:-1]
        out_list.append(unquoted)
    return out_list


def parse_clade_defs(clade_defs_fp):
    if not clade_defs_fp:
        return None
    clades = {}
    with open(clade_defs_fp, "r") as inp:
        for line in inp:
            line_stripped = line.strip()
            try:
                name, defn = line_stripped.split("\t")
            except:
                m = f"Expecting 2 tab separated columns found {repr(line)} in {clade_defs_fp}"
                raise RuntimeError(m)
            if name in clades:
                raise RuntimeError(f"Clade name {name} repeated.")
            clades[name] = _parse_cdef_str(defn)
    return clades
