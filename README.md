# SEVA
Code for analyses in [Splice Expression Variation Analysis (SEVA) for Inter-tumor Heterogeneity of Gene Isoform Usage in Cancer](https://doi.org/10.1101/091637)

**Authors:** Bahman Afsari, Theresa Guo, Michael Considine, Liliana Florea, Dylan Kelley, Emily Flam, Patrick K. Ha, Donald Geman, Michael F. Ochs, Joseph A. Califano, Daria A. Gaykalova, Alexander V. Favorov, Elana J. Fertig

## Abstract

**Motivation:** Alternative splicing events (ASE) cause expression of a variable repertoire of potential protein products that are critical to carcinogenesis. Current methods to detect ASEs in tumor samples compare mean expression of gene isoforms relative to that of normal samples. However, these comparisons may not account for heterogeneous gene isoform usage that is common in tumors.

**Results:** Therefore, we introduce Splice Expression Variability Analysis (SEVA) to detect differential splice variation, which accounts for tumor heterogeneity. This algorithm compares the degree of variability of junction expression profiles within a population of normal samples relative to that in tumor samples using a rank-based multivariate statistic that models the biological structure of ASEs. Simulated data show that SEVA is more robust to tumor heterogeneity and its candidates are more independent of differential expression than EBSeq, DiffSplice, and rMATS. SEVA analysis of head and neck tumors identified differential gene isoform usage robust in cross-study validation. The algorithm observes substantial differences in gene isoform usage between head and neck tumor subtypes, with greater inter-tumor heterogeneity in HPV-negative tumors with alterations to genes that regulate RNA splice machinery. Thus, SEVA is well suited for differential ASE analysis and assessment of ASE inter-tumor heterogeneity in RNA-seq data from primary tumor samples.

**Availability:** SEVA is implemented in the R/Bioconductor package [GSReg](https://bioconductor.org/packages/release/bioc/html/GSReg.html).

## Contact

* [Bahman Afsari](bahman@jhu.edu)
* [Daria Gaykalova](dgaykal1@jhmi.edu)
* [Alexander Favorov](favorov@sensi.org)
* [Elana Fertig](ejfertig@jhmi.edu)
