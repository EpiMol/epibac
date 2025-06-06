# Changelog

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
