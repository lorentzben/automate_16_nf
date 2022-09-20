#!/usr/bin/bash

eval "$(conda shell.bash hook)"

conda activate python2

item=$(sed -n '1p' current.txt)

export2graphlan.py -i $item\-otu-table-mod.biom -t tree.nwk -a annot.txt --discard_otus --most_abundant 40  --annotations 2,3,4,5 --external_annotations 6 --def_na 0 --internal_levels --title "OTU Tree of "$item" by GraPhlAn"

graphlan_annotate.py --annot annot.txt tree.nwk annotated_tree.txt

graphlan.py annotated_tree.txt $item\_image_graph.png --format png --size 10 --dpi 100

graphlan.py annotated_tree.txt $item\_image__pdf_graph.png --format png  --dpi 72

conda deactivate