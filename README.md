R script (civic.R) to tidy data on prognostic/diagnostic/drug sensitivity biomarkers for personalized oncology

Brief overview:
* Variants are retrieved through the API of [Clinical Interpretation of Variants in Cancer](https://genome.github.io/civic-api-docs/) - data fetched 20160917
* HTML links are created for sources of evidence (i.e. PubMed IDs) that underlie individual markers
* The following biomarker types are currently kept:
  - Mutations (missense/stop_gained etc.)
  - Insertion/deletions, frameshifts
  - Gene codons
  - Gene exon/domain mutations
  - Gene copy number aberrations (loss, gain)
  - Expression
* Biomarkers that lack a mapping to the genome (GRCh37) are omitted
* Two output files are created:
  - __civic.biomarkers.tsv__ : Tab-separated file with data (gene, drugs, association etc.) on individual biomarkers
  - __civic.bed__ : BED file with all chromosomal regions (and identifiers) of biomarkers in _civic.biomarkers.tsv_
