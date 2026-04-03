Bundled `inst/extdata` provenance
=================================

This directory documents how the files shipped in `inst/extdata/` were
generated and where the underlying source material came from.

Files currently bundled
-----------------------

1. `inst/extdata/foldEnergies/human/humanDB_UTR5_foldEnergy.txt.gz`
2. `inst/extdata/foldEnergies/human/humanDB_CDS_foldEnergy.txt.gz`
3. `inst/extdata/foldEnergies/human/humanDB_UTR3_foldEnergy.txt.gz`
4. `inst/extdata/foldEnergies/mouse/mouseDB_UTR5_foldEnergy.txt.gz`
5. `inst/extdata/foldEnergies/mouse/mouseDB_CDS_foldEnergy.txt.gz`
6. `inst/extdata/foldEnergies/mouse/mouseDB_UTR3_foldEnergy.txt.gz`
7. `inst/extdata/indexes/human/IndexesHuman.txt`
8. `inst/extdata/indexes/mouse/IndexesMouse.txt`

Overview of source data
-----------------------

The bundled files are derived summary tables produced from reference
transcript/coding-sequence annotations used by `postNet`.

- Human and mouse transcript annotations used by `postNet` are based on NCBI
  RefSeq releases described in `?postNetStart`.
- Codon-related summaries use coding sequences corresponding to the human and
  mouse reference annotations used by the package and/or NCBI CCDS resources as
  described in `?codonUsage`.
- Folding-energy tables are derived from the same packaged reference transcript
  regions and were calculated externally with the mfold/UNAFold algorithm, as
  described in `?foldingEnergyAnalysis`.

How the `foldEnergies` files were generated
-------------------------------------------

Each gzipped table contains three columns:

- `id`: transcript identifier
- `fold_energy`: predicted Gibbs free energy
- `length`: sequence length for that transcript region

High-level procedure:

1. Start from the packaged human or mouse transcript reference annotation used
   by `postNet` (`UTR5`, `CDS`, and `UTR3` regions).
2. For each transcript region sequence, run RNA folding with the mfold/UNAFold
   algorithm.
3. Extract the reported folding energy for the selected transcript region.
4. Store one row per transcript as a tab-delimited text file.
5. Compress the final table with gzip.

Pseudo-code:

```text
for species in {human, mouse}:
    for region in {UTR5, CDS, UTR3}:
        seqs <- reference annotation for species/region
        for each transcript sequence:
            run mfold/UNAFold on the sequence
            record transcript id, predicted free energy, sequence length
        write tab-delimited table
        gzip table
```

How the `indexes` files were generated
--------------------------------------

Each index table contains one row per gene and the columns:

- `external_gene_name`
- `CAI`
- `CBI`
- `Fop`
- `L_aa`
- `tAI`

High-level procedure:

1. Obtain coding sequences for the supported species from the reference coding
   sequence set used by `postNet`.
2. Translate or otherwise summarize each coding sequence at the gene level.
3. Compute the codon-usage metrics used by the package:
   `CAI`, `CBI`, `Fop`, `L_aa`, and `tAI`.
4. Write one tab-delimited summary table per species.

Pseudo-code:

```text
for species in {human, mouse}:
    cds <- reference coding sequences for species
    group sequences by gene
    for each gene:
        compute CAI
        compute CBI
        compute Fop
        compute protein length / L_aa
        compute tAI
    write tab-delimited summary table
```

Source and licensing notes
--------------------------

- NCBI RefSeq: source annotations are from the NCBI Reference Sequence
  database. See <https://www.ncbi.nlm.nih.gov/refseq/>.
- NCBI CCDS: coding-sequence resources referenced by `postNet` are from the
  Consensus CDS project. See
  <https://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi>.
- mfold/UNAFold: folding-energy values were generated with the mfold/UNAFold
  software described in `?foldingEnergyAnalysis`. See
  <https://www.unafold.org/>.

Licensing/status of bundled material:

- The files in `inst/extdata/` are derived summary tables created for use by
  this package rather than verbatim redistribution of the complete upstream
  source databases.
- NCBI resources are public reference resources provided by the U.S.
  National Library of Medicine / NCBI; users should consult the upstream NCBI
  terms and citation guidance for reuse of the original source data.
- mfold/UNAFold is external software; users should consult the upstream site
  for software terms, citation requirements, and any usage restrictions.

Maintenance note
----------------

If any bundled file in `inst/extdata/` is replaced or regenerated, this file
should be updated to record the exact upstream release, date, and command or
script used for regeneration.
