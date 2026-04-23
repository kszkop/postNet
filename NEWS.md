--------------------------

# postNet ver. 0.99.9
## 2026-04-23

    - Fixed feature integration network plotting when feature-feature associations
      are missing or undefined, including cases where no associations pass the
      `covarFilt` threshold or correlations are `NA`. The network plot now renders
      the selected feature nodes without unintended default igraph edges, and
      grouped-layout helper edges are used only for layout calculation rather than
      being retained in the plotted network.
      
    - Fixed `uorfAnalysis()` position output when `onlyUTR5 = TRUE` and
      `unitOut = "position"`, allowing genes with multiple uORFs to return
      start and end positions without vector-length errors.

--------------------------

# postNet ver. 0.99.0
## 2026-01-26

    - Initial submission to Bioconductor
