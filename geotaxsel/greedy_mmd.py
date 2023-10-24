#! /usr/bin/env python3
from dendropy.calculate.phylogeneticdistance import PhylogeneticDistanceMatrix
from geopy.distance import geodesic

from geotaxsel import debug, info


def calc_dist(loc_1, loc_2):
    if loc_1 is loc_2:
        return 0.0
    return geodesic(loc_1.coords, loc_2.coords).km


def sel_most_geo_div_taxon(label_ind_pairs, loc_list, sp_by_name):
    md = None
    md_loc = None
    md_sp_ind = None
    for label, ind in label_ind_pairs:
        sp = sp_by_name[label]
        for l1 in sp.locations:
            sum_sq_dist = 0.0
            for l2 in loc_list:
                try:
                    sum_sq_dist += calc_dist(l1, l2)
                except:
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


def tip_to_root_dist(nd, root):
    t = 0.0
    while nd is not root:
        t += nd.edge_length
        nd = nd.parent_node
        assert nd is not None
    return t


def ultrametric_greedy_mmd(tree, num_taxa, sp_by_name):
    tree.calc_node_ages(ultrametricity_precision=5e-5)
    # for nd in tree.ageorder_node_iter(descending=True):
    #     if nd.is_leaf():
    #         print(f"{nd.taxon.label} -> root = ", tip_to_root_dist(nd, tree.seed_node))
    #     else:
    #         print(
    #             f"internal node at age {nd.age} and tip to root =",
    #             tip_to_root_dist(nd, tree.seed_node),
    #         )
    last_added = {tree.seed_node}
    chosen_ancs = set(last_added)
    for nd in tree.ageorder_node_iter(descending=True):
        if len(chosen_ancs) >= num_taxa:
            break
        if nd.is_leaf():
            raise NotImplementedError(
                "num_taxa is greater than the number of internal nodes. Not implemented yet."
            )
        assert nd in chosen_ancs
        chosen_ancs.remove(nd)
        last_added = set(nd.child_nodes())
        chosen_ancs.update(last_added)
    if len(chosen_ancs) != num_taxa:
        raise NotImplementedError(
            "Polytomy caused num_taxa to be exceeded need to check last_added and remove some..."
        )
    return chosen_ancs


def subtree_name(nd, mrca_notation=True):
    if hasattr(nd, "clade_names"):
        return nd.clade_names[-1]
    if nd.is_leaf():
        return nd.taxon.label
    children = nd.child_nodes()
    assert len(children) > 1
    child_names = [subtree_name(i, mrca_notation=False) for i in children]
    child_names_w_len = [(len(i), i) for i in child_names]
    shortest, sec_shortest = child_names_w_len[0][1], child_names_w_len[1][1]
    if mrca_notation:
        return f"MRCA({shortest} + {sec_shortest})"
    else:
        return shortest


def anc_name(nd):
    par = nd.parent_node
    if par is None:
        return ""
    if hasattr(par, "clade_names"):
        return par.clade_names[0]
    return anc_name(par)


def greedy_mmd(tree, num_taxa, sp_by_name):
    """Selects `num_taxa` from the tips of `tree`.

    `sp_by_name` should be a dict mapping a name to Species object.
    Every tip label in the tree must be in sp_by_name
    """
    # not efficiently updating min dist yet
    taxa_list = [i.taxon for i in tree.leaf_nodes()]
    taxa_label_list = [i.label for i in taxa_list]
    if num_taxa == len(taxa_label_list):
        return taxa_label_list
    if num_taxa > len(taxa_label_list):
        raise ValueError("num_taxa exceeds the number of taxa in the tree")
    debug(f"Calculating patristic distance matrix")
    pdmc = PhylogeneticDistanceMatrix.from_tree(tree)
    # Making a shallow copy of the distance matrix here. just being lazy...
    # While copying, find the pair of taxa with the largest patristic distance.
    # TODO: this does not consider ties in the largest distance
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
            ndist = pdmc.patristic_distance(row_tax, col_tax)
            if max_dist is None or ndist > max_dist:
                max_dist = ndist
                max_dist_inds = (row_ind, col_ind)
            row.append(ndist)
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
            assert inc_dist
            curr_min_dist = min(inc_dist)
            if curr_min_dist > max_min_dist:
                max_min_dist = curr_min_dist
                mmd_ind_set = set()
                mmd_ind_set.add(row_ind)
            elif abs(curr_min_dist - max_min_dist) < TOL:
                mmd_ind_set.add(row_ind)
        tied_tax = [(taxa_label_list[i], i) for i in mmd_ind_set]
        mmd_ind, sel_loc = sel_most_geo_div_taxon(tied_tax, curr_locs, sp_by_name)
        ntl = taxa_label_list[mmd_ind]
        debug(
            f'    taxon {1 + len(sel_inds)} = "{ntl}" with MD = {max_min_dist} set = {mmd_ind_set}'
        )
        sel_tax_labels.append(ntl)
        sel_inds.add(mmd_ind)
        curr_locs.append(sel_loc)
    debug("Taxa selected, cleaning up...")
    return sel_tax_labels


def output_chosen_anc(tree, cut_branches_fp, chosen_ancs):
    if not cut_branches_fp:
        return
    single_mode = True
    if tree is None:
        single_mode = False
        header = "num-leaves\tnum-times-selected\tnames\n"
    else:
        header = (
            "dist-to-root\tdist-to-tips\tnum-leaves\tname-or-mrca\tancestral-clade\n"
        )

    if single_mode:
        with open(cut_branches_fp, "w") as outp:
            outp.write(header)
            prod = 1
            for nd in chosen_ancs:
                t2rd = tip_to_root_dist(nd, tree.seed_node)
                leaves_below = list(nd.leaf_nodes())
                nl = len(leaves_below)
                prod *= nl
                nm = subtree_name(nd)
                ancnm = anc_name(nd)
                outp.write(f"{t2rd}\t{nd.age}\t{nl}\t{nm}\t{ancnm}\n")
            prod_str = str(prod)
            pow_10 = len(prod_str) - 1
            if pow_10 > 9:
                cut_prod_str = prod_str[:7]
                dec_str = f"{cut_prod_str[0]}.{cut_prod_str[1:]}"
                parens_str = f" (about {dec_str}E{pow_10})"
            else:
                parens_str = ""
            msg = f"{prod_str}{parens_str} total choices of {len(chosen_ancs)} that maximize PD."
            info(msg)
    else:
        with open(cut_branches_fp, "w") as outp:
            outp.write(header)
            for k, v in chosen_ancs.items():
                lk = list(k)
                lk.sort()
                lf = ", ".join(lk)
                outp.write(f"{len(lk)}\t{v}\t{lf}\n")
