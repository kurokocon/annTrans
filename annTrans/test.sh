#!/bin/bash

./anno_transfer --target fixed.gff3 --source tri6.gene_structures_post_PASA_updates.65674.gff3 --output test --genome-sequence genome_S288C_R64.fsa --target-overlap 0.9 --source-overlap 0.9
