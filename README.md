# taxon-selection
Code for selecting taxa that are geographically and phylogenetically 
dissimilar.

Installation and venv setup:

    python -mvenv env
    source env/bin/activate
    pip install DendroPy
    pip install geopy
    python setup.py

Data setup:

    mkdir cruft
    cd cruft
    wget https://www.mammaldiversity.org/assets/data/MDD.zip
    unzip MDD.zip
    wget https://github.com/n8upham/MamPhy_v1/raw/master/_DATA/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_MCC_v2_target.tre
    cd ..
    mkdir data

That tree was linked from https://vertlife.org/data/mammals/

Then copy `outfile_rmDups_filtered_4Holder_v4.csv` into `data`

Every session:

    source env/bin/activate

Creating clades file for updating tip names and naming internals:

    python taxonomy-to-clades.py \
        cruft/MDD/MDD_v1.11_6649species.csv \
        cruft/cmw2mdd.tsv > cruft/clades.tsv 

Running taxon selection (currently will crash before final collection with an "early" exit):

    python taxselect.py \
        --centroid-file=data/outfile_rmDups_filtered_4Holder_v4.csv \
        --tree-file=cruft/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_MCC_v2_target.tre \
        --num-to-select=N \
        --clade-defs-file=cruft/clades.tsv \
        --cut-branches-file=cruft/cut_branches.csv \
        --name-updating-file=cruft/cmw2mdd.tsv


where `N` is the number of taxa to choose.

Has been tested with DendroPy-4.6.1 and geopy-2.4.0 and Python 3.10.12 on Ubuntu

## Info in the current error stream
  * "Skipping taxon "X" due to NA in centroid." indicates that `X` was not found in the centroid-file
  * "tip names updated to new taxonomy" indicates use of the "CMW_sciName" field of the MDD taxonomy to update taxon 
     labels in the tree from `CMW_sciName` to MDD's `sciName` to improve matching.
  * "Will prune X: not matching expected form of a species name" refers to tree tips removed for not looking like
     species names
  * "Will prune X: not found in data/outfile_rmDups_filtered_4Holder_v4.csv" refers to tips in the tree
     not found in the centroid file. So these are pruned. *Pruning may not be the right choice*
  * "Will prune X: not found in any clade in clade definitions" tips not found in MDD. that make it tough to
    use clade names in output. **We probably *DON'T* want to prune these in the final analysis**
  * "X taxa in clade definitions, but not in the tree:" refers to species in MDD but not the tree.
  * then comes a list of taxa in MDD that are not monophyletic in the tree
  * then a list of taxa in MDD that are monophyletic.