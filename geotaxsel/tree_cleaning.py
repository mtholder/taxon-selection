#! /usr/bin/env python3
import re
import itertools
from .logs import info


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


def prune_taxa_without_sp_data(
    tree,
    sp_w_data,
    upham_to_iucn=None,
    name_mapping_fp="",
    centroid_fp="",
    clades=None,
    new_names_for_leaves=None,
):
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
    sys.exit("early")
