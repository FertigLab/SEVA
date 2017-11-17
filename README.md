# SEVA
Code for analyses in [Splice Expression Variation Analysis (SEVA) for Inter-tumor Heterogeneity of Gene Isoform Usage in Cancer](https://doi.org/10.1101/091637)

**Authors:** Bahman Afsari, Theresa Guo, Michael Considine, Liliana Florea, Luciane T. Kagohara, Dylan Kelley, Emily Flam, Patrick K. Ha, Donald Geman, Michael F. Ochs, Genevieve L. Stein-O'Brien, Joseph A. Califano, Daria A. Gaykalova, Alexander V. Favorov, Elana J. Fertig

## Abstract

**Motivation:** Current bioinformatics methods to detect changes in gene isoform usage in distinct phenotypes compare the relative expected isoform usage in phenotypes. These statistics accurately model the difference in isoform usage in normal tissues, which have stable regulation of gene splicing. Pathological conditions, such as cancer, can have broken regulation of splicing that induces changes in the heterogeneity of the expression of splice variants between samples rather than changes in relative expected expression. Inferring events with such differential heterogeneity in gene isoform usage requires new statistical approaches.

**Results:** We introduce Splice Expression Variability Analysis (SEVA) to model increased heterogeneity of splice variant usage between conditions (e.g., tumor and normal samples). SEVA uses a rank-based multivariate statistic that compares the variability of junction expression profiles within one condition to the variability within another.  Simulated data show that SEVA is unique in modeling heterogeneity of gene isoform usage, and benchmark SEVA's performance against the differential gene isoform usage algorithms EBSeq, DiffSplice, and rMATS that do not explicitly model this heterogeneity. We confirm the accuracy of SEVA in identifying known splice variants in head and neck cancer and perform cross-validation of novel splice variants learned in both sample cohorts. A novel comparison of splice variant heterogeneity between subtypes of head and neck cancer demonstrated unanticipated similarity between the heterogeneity of gene isoform usage in HPV-positive and HPV-negative subtypes and anticipated increased heterogeneity among HPV-negative samples with mutations in genes that regulate the splice variant machinery.

**Conclusion:** These results show that SEVA accurately models differential heterogeneity of gene isoform usage from RNA-seq data.

**Availability:** SEVA is implemented in the R/Bioconductor package [GSReg](https://bioconductor.org/packages/release/bioc/html/GSReg.html).

## Contact

* [Bahman Afsari](bahman@jhu.edu)
* [Alexander Favorov](favorov@sensi.org)
* [Elana Fertig](ejfertig@jhmi.edu)
