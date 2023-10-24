#! /usr/bin/env python3
import re
import itertools
from .logs import info
from .taxonomy import Ranks


def tips_from_clades(clades):
    leaf_label_set = set()
    for name, cd in clades.items():
        for label in itertools.chain(cd.must, cd.might):
            if label not in clades:
                leaf_label_set.add(label)
    info(f"{len(leaf_label_set)} labels in clades")
    return leaf_label_set


def find_tree_leaves_not_in_clades(clades, tree):
    clade_tips = tips_from_clades(clades)
    to_del = []
    for leaf in tree.leaf_nodes():
        if not leaf.taxon.label in clade_tips:
            print(leaf.taxon.label)
            to_del.append(leaf)
    return to_del


def alert_pruning(msg, to_prune):
    if not to_prune:
        return
    info(f"Will prune {len(to_prune)}: {msg}")
    for el in to_prune:
        info(f'  "{el}"')


def _fill_must_taxa(cdef, missing_from_tree, tip_label_2_nd, clade_name_to_def, strict):
    reprocess = False
    for label in cdef.must:
        if label in missing_from_tree:
            continue
        nd = tip_label_2_nd.get(label)
        if nd is None:
            des_cdef = clade_name_to_def.get(label)
            if des_cdef is None:
                if len(label.split()) == 2:
                    missing_from_tree.add(label)
                elif strict:
                    raise RuntimeError("Higher taxon {label} not found")
                else:
                    reprocess = True
            else:
                cdef.must_taxa.update(des_cdef.must_taxa)
        else:
            cdef.must_taxa.add(nd.taxon)
    return reprocess


def _fill_might_taxa(
    cdef, missing_from_tree, tip_label_2_nd, clade_name_to_def, strict
):
    reprocess = False
    for label in cdef.might:
        if label in missing_from_tree:
            continue
        nd = tip_label_2_nd.get(label)
        if nd is None:
            des_cdef = clade_name_to_def.get(label)
            if des_cdef is None:
                if len(label.split()) == 2:
                    missing_from_tree.add(label)
                elif strict:
                    raise RuntimeError(f"Higher taxon {label} not found")
                else:
                    reprocess = True
            else:
                cdef.might_taxa.update(des_cdef.must_taxa)
                cdef.might_taxa.update(des_cdef.might_taxa)
        else:
            cdef.might_taxa.add(nd.taxon)
    return reprocess


def link_cdef_to_taxa(cdefs_by_rank, tip_label_2_nd):
    missing_from_tree = set()
    clade_name_to_def = {}
    for cdef_list in cdefs_by_rank:
        reprocess_list = []
        for cdef in cdef_list:
            cdef.must_taxa = set()
            if _fill_must_taxa(
                cdef, missing_from_tree, tip_label_2_nd, clade_name_to_def, strict=False
            ):
                reprocess_list.append(cdef)
            cdef.might_taxa = set()
            if _fill_might_taxa(
                cdef, missing_from_tree, tip_label_2_nd, clade_name_to_def, strict=False
            ):
                reprocess_list.append(cdef)
            clade_name_to_def[cdef.name] = cdef
        for cdef in reprocess_list:
            _fill_must_taxa(
                cdef, missing_from_tree, tip_label_2_nd, clade_name_to_def, strict=True
            )
            _fill_might_taxa(
                cdef, missing_from_tree, tip_label_2_nd, clade_name_to_def, strict=True
            )
    return missing_from_tree


def _eval_one_clade(tree, cdef, clades, tip_label_2_nd):
    if not cdef.must_taxa:
        return False, None
    one_taxon = next(iter(cdef.must_taxa))
    one_leaf = tip_label_2_nd[one_taxon.label]
    curr_nd = one_leaf
    while True:
        nls = set(curr_nd.bipartition.leafset_taxa(tree.taxon_namespace))
        nld = nls.difference(cdef.must_taxa)
        if not nld.issubset(cdef.might_taxa):
            nldt = [i.label for i in nld]
            if cdef.must_taxa.issubset(nls):
                info(f"Clade {cdef.name} not found due to intrusion of {nldt}")
            else:
                lsd = nls.difference(cdef.must_taxa)
                lsdt = [i.label for i in lsd]
                lsi = nls.intersection(cdef.must_taxa)
                lsit = [i.label for i in lsi]
                info(
                    f"Clade {cdef.name} not found due to intrusion of {nldt} as closer to {lsit} than the other members({lsdt}) are."
                )
            return False, curr_nd
        if cdef.must_taxa.issubset(nls):
            return True, curr_nd
        curr_nd = curr_nd.parent_node
        assert curr_nd is not None


def _link_existing_clades(tree, cdefs_by_rank, clades, tip_label_2_nd):
    found, not_found = {}, {}
    for cdef_list in cdefs_by_rank:
        for cdef in cdef_list:
            is_in_tree, nd = _eval_one_clade(tree, cdef, clades, tip_label_2_nd)
            if is_in_tree:
                found[cdef.name] = nd
                cdef.node = nd
                if hasattr(nd, "clade_names"):
                    nd.clade_names.append(cdef.name)
                else:
                    nd.clade_names = [cdef.name]
            else:
                not_found[cdef.name] = nd
                cdef.conflicting = nd
    return found, not_found


def label_internals(tree, clades):
    tip_label_2_nd = {}
    for leaf in tree.leaf_node_iter():
        assert leaf.taxon.label not in tip_label_2_nd
        tip_label_2_nd[leaf.taxon.label] = leaf
    max_rank_v = max([i.value for i in Ranks])
    cdefs_by_rank = [[] for i in range(1 + max_rank_v)]
    for name, cdef in clades.items():
        cdef.name = name
        assert cdef.rank  # May not always hold, but should for mammal taxonomy from MDD
        cdefs_by_rank[cdef.rank.value].append(cdef)
    cdefs_by_rank.reverse()
    missing_from_tree = link_cdef_to_taxa(cdefs_by_rank, tip_label_2_nd)
    if missing_from_tree:
        info(
            f"{len(missing_from_tree)} taxa in clade definitions, but not in the tree:"
        )
        for el in missing_from_tree:
            info(f'  "{el}"')
    found, not_found = _link_existing_clades(
        tree, cdefs_by_rank, clades, tip_label_2_nd
    )
    info(f"{(len(not_found))} clades broken.")
    info(f"{(len(found))} clades found:")
    for c in found:
        info(f"  {c} found")


def prune_taxa_without_sp_data(
    tree,
    sp_w_data,
    upham_to_iucn=None,
    name_mapping_fp="",
    centroid_fp="",
    clades=None,
    new_names_for_leaves=None,
):
    info(
        f"prune_taxa_without_sp_data(tree, sp_w_data, upham_to_iucn={upham_to_iucn}...)"
    )
    clade_tips = tips_from_clades(clades) if clades else set()
    if new_names_for_leaves is None:
        new_names_for_leaves = {}
    taxa_list = [i for i in tree.taxon_namespace]
    to_prune = []
    sp_pat = re.compile("^([A-Z][a-z]+ +[-a-z0-9]+) [A-Z][A-Za-z]+ [A-Z]+$")
    bad_names = []
    no_geo = []
    null_name_mapped = []
    not_in_clades = []
    remapped = []
    final_name_set = set()
    for n, i in enumerate(taxa_list):
        m = sp_pat.match(i.label)
        if not m:
            to_prune.append(i)
            bad_names.append(i.label)
            continue
        sp_name = m.group(1)
        if upham_to_iucn is not None:
            next_name = upham_to_iucn.get(sp_name)
        else:
            next_name = sp_name
        if next_name is None:
            to_prune.append(i)
            null_name_mapped.append(sp_name)
            continue
        final_name = new_names_for_leaves.get(next_name, next_name)
        if final_name != next_name:
            remapped.append((next_name, final_name))
        if final_name not in sp_w_data:
            to_prune.append(i)
            no_geo.append(final_name)
        elif clades and final_name not in clade_tips:
            to_prune.append(i)
            not_in_clades.append(final_name)
        else:
            final_name_set.add(final_name)
            i.label = final_name
    info(f"{len(remapped)} tip names updated to new taxonomy:")
    for from_n, to_n in remapped:
        info(f"  {from_n} --> {to_n}")
    alert_pruning(f"not matching expected form of a species name", bad_names)
    alert_pruning(f"not found in {name_mapping_fp}", null_name_mapped)
    alert_pruning(f"not found in {centroid_fp}", no_geo)
    alert_pruning(f"not found in any clade in clade definitions", not_in_clades)
    tree.prune_taxa(to_prune)
    if clades:
        tree.encode_bipartitions()
        label_internals(tree, clades)
    centroids_but_no_tips = set()
    for sp_name in sp_w_data:
        if sp_name not in final_name_set:
            centroids_but_no_tips.add(sp_name)
    info(f"{len(centroids_but_no_tips)} species in centroid file but not in the tree.")
    for sp_name in centroids_but_no_tips:
        info(f"  {sp_name}")
