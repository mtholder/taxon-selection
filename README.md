# taxon-selection
Code for selecting taxa that are geographically and phylogenetically 
dissimilar.

    python -mvenv env
    source env/bin/activate
    pip install DendroPy
    pip install geopy
    python taxselect.py \
        --centroid-file=path-to-centroids.csv \
        --tree-file=ultrametric.tre \
        --num-to-select=N

where `N` is the number of taxa to choose.

Has been tested with DendroPy-4.6.1 and geopy-2.4.0 and Python 3.10.12 on Ubuntu