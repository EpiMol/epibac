# Changelog

## [1.2.3] - 2026-06-03

### Fixed
- `mob_suite=3.1.9` no longer breaks the creation of the `epibac_amr_annotation`
  env because we pin `biopython=1.81` (required for `gc_fraction`, added in
  Biopython ≥ 1.80).
- `mob_suite=3.1.9` no longer fails at runtime with
  `TypeError: NDFrame.to_csv() got an unexpected keyword argument 'line_terminator'`
  because we pin `pandas=1.5.3` in the 3 envs that use mob_suite
  (`epibac_amr_annotation`, `epibac_mge_mob_suite`, `epibac_mobile_elements`).
  The kwarg was removed in pandas 2.0.
- `epibac.py` now launches Snakemake with `PYTHONNOUSERSITE=1`. This isolates the
  conda envs from the operator's `~/.local/lib/pythonX.Y/site-packages/`, preventing
  a `pip install --user` (e.g. `pandas-2.0.1`) from silently overriding the
  packages pinned in each rule's yml.

### Documentation
- README: new Section 2 "Quickstart — routine use" for hospital staff,
  with the 4 commands for each run, run naming, and samplesheet format.

## [1.2.2] - 2026-05-06

> Backfilled retroactively — the tag was created without a CHANGELOG entry.

### Added
- RGI (CARD) AMR detection: new conda env (`epibac_rgi.yml`), `setup_rgi_database`
  rule, `epibac_rgi` analysis rule, config params/resources, and inclusion in
  pipeline outputs via `get_all_inputs()`.

### Changed
- Kraken2 standard DB URL updated to the 16 GB / 2025-07-14 build.
- `kraken2` pinned to `2.1.6` in the QC env.
- `epibac_summary_gestlab.py`: plasmid summary rewritten to read per-sample
  MOB-suite `mobtyper_results.txt` and produce detailed per-plasmid lines;
  added companion `*_Plasmids_Detail.xlsx` report.

### Fixed
- `validate_samples_file.py`: detect and strip UTF-8 BOM (Excel-on-Windows),
  smarter `;` vs `,` separator detection, and duplicate sample-ID handling.
  Always save the validated CSV as UTF-8 without BOM.
- `epibac.py` / `common.smk`: write validated CSV without BOM and auto-detect
  the separator when reloading it downstream.

## [1.2.1] - 2025-06-10

> Backfilled retroactively — the tag was created without a CHANGELOG entry.

### Fixed
- Platon results capture: dynamic variables instead of hardcoded paths for the
  Platon database, and corrected JSON file pattern matching in
  `unified_plasmid_analysis.py` so all `.json` files are found.
- Numeric sample IDs (e.g. `234512`, `234518`, `234656`) now work end-to-end:
  sample-ID extraction logic improved, and int-vs-string comparison fixed in
  the GESTLAB report generator.
- MOB-suite and Platon integration no longer fail for numeric-only sample IDs.

### Documentation
- README cleanup: removed malformed HTML and "Copiar" buttons left over from
  browser copy-paste, all links converted to proper Markdown.
- Added the complete list of 21 supported bacterial species under
  `ESPECIE_SECUENCIA`.
- Practical example: show all 3 samples (`234512`, `234518`, `234656`) and
  point to the Zenodo dataset `250425_GRAL001`.

## [1.2.0] - 2025-06-06

### Enhanced MLST Analysis
- Integrated custom MLST database for all supported species, replacing pubMLST
- Built custom BLAST database for MLST detection using assembled allele files
- Added MLST schemes: Citrobacter spp., E. coli Achtman, K. oxytoca
- Reviewed scheme-exclusion logic to avoid unintended filters
- Fixed species-name conversion based on MLST results for GestLab export
- Corrected extensions and formatting of MLST schemas

### Nanopore & Dorado Integration  
- Added per-sample Dorado model specification for hybrid assemblies
- Support for flexible basecalling models (fast/hac/sup) and versions (4.3/5.0)
- Enhanced Nanopore sample processing pipeline
- Updated sample validation schema with 'dorado_model' column
- Improved hybrid assembly workflows

### Mobile Genetic Elements Analysis
- Added Platon database download rule for plasmid analysis
- Implemented Snakemake rules for Platon and MOB-suite
- Custom parsers for Platon JSON and MOB-suite TSV
- Integrated plasmid data with AMR and virulence outputs
- Prepared transposon database for Abricate
- Added conda environments for MGE analysis

### GestLab Report Enhancements
- Enhanced output format with improved consistency
- Collapsed plasmid results into readable PLASMIDOS_WGS column
- Comprehensive plasmid summaries (count, size, rep types, mobility)
- Fixed Excel export compatibility with pandas NA values
- Improved error handling in report generation

### Technical Improvements
- Enhanced validation system and Snakemake command handling
- Improved setup and configuration management
- Better database management with automated procedures
- Enhanced logging and debugging capabilities
- Streamlined workflow integration

### Database Updates
- Custom MLST database for all supported species
- Transposon database with resistance-associated elements
