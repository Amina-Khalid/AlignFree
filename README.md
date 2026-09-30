# AlignFree

A GUI-based tool for alignment-free sequence comparison and phylogenetic reconstruction.

## What it does

AlignFree compares biological sequences (e.g. genes) without requiring sequence alignment.
It represents each sequence as a k-mer (short-word) frequency profile, computes Euclidean
distance between sequences, and reconstructs a phylogenetic tree (UPGMA or Neighbor-Joining)
from the resulting distance matrix.

## Current features

- Pairwise sequence comparison (two named sequences, k-mer length, Euclidean distance)
- Multiple sequence comparison (distance matrix across any number of sequences)
- UPGMA / Neighbor-Joining phylogenetic tree construction and visualization

## Setup

```bash
pip install -r requirements.txt
python AlignFree.py
```

## Status

Early/initial version. See commit history for ongoing development.

## License

All rights reserved. This code is made publicly viewable for academic
review and demonstration purposes only. See the [LICENSE](LICENSE) file
for details. Reuse, redistribution, or publication of this code or any
derivative of it is not permitted without prior written permission from
the author.

