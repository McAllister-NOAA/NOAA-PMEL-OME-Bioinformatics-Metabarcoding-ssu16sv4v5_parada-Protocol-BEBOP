---
# MIOP terms
methodology category: Omics Analysis
project: NOAA Pacific Marine Environmental Laboratory Ocean Molecular Ecology Group protocols
purpose: 'taxonomic diversity assessment by targeted gene survey [OBI:0001960]'
analyses: 'bioinformatics analysis [MIOP:0000004] | amplicon sequencing assay [OBI:0002767]'
geographic location: 
   options: 'North East Pacific Ocean [GAZ:00013765] | Bering Sea [GAZ:00008990] | Arctic Ocean [GAZ:00000323]'
broad-scale environmental context: 'marine biome [ENVO:00000447] | oceanic epipelagic zone biome [ENVO:01000035]'
local environmental context: 
   options: 'marine photic zone [ENVO:00000209] | marine aphotic zone [ENVO:00000210] | marine benthic biome [ENVO:01000024]'
environmental medium: 'ocean water [ENVO:00002149] | sea water [ENVO:00002149]'
target: 'Bacterial 16S RNA [NCIT:C105370]'
creator: Samantha Setta, Sean McAllister, Zachary Gold
materials required: high-performance computing resources
skills required: basic bash, R
time required: variable
personnel required: 1
language: en
issued: 2026-01-27
audience: scientists
publisher: 'NOAA Pacific Marine Environmental Laboratory Ocean Molecular Ecology Group | University of Washington Cooperative Institute for Climate, Ocean, & Ecosystem Studies'
hasVersion: 1
license: 'CC0 1.0 Universal'
maturity level: demonstrated

# FAIRe terms
sop_bioinformatics: TBD_this_DOI (paste link when published)
checkls_ver: TBD_bioinformatics_template
assay_name: 'ssu16sv4v5_parada | ssu16sv4v5_parada_OSUmod'
pcr_primer_forward: GTGYCAGCMGCCGCGGTAA
pcr_primer_reverse: CCGYCAATTYMTTTRAGTTT
trim_method: 'Cutadapt, primer trimming | DADA2, filterAndTrim (quality and length trimming)'
trim_param:
   - Cutadapt:
      default: '-a "GTGYCAGCMGCCGCGGTAA;required...TTACCGCGGCKGCTGRCAC;optional", -A "CCGYCAATTYMTTTRAGTTT;required...AAACTYAAAKRAATTGRCGG;optional", --discard-untrimmed, -m 1'
   - DADA2:
      default: 'trunQ = {dada_trunQ}, trimRight = {dada_trimRight}, trimLeft = {dada_trimLeft}'
      source_file: REVAMP_config
      source_term: 'dada_trunQ | dada_trimRight | dada_trimLeft'
demux_tool: 'pheniqs v2.1.0 | Cutadapt v3.4'
demux_max_mismatch: '0 | 0.3'
merge_tool: 'DADA2, mergePairs'
merge_min_overlap: 20
min_len_cutoff:
   default: 100
   source_file: REVAMP_config
   source_term: dada_minlength
min_len_tool: DADA2
error_rate_tool: DADA2
error_rate_type: expected error rate
error_rate_cutoff:
   default: 2
   source_file: REVAMP_config
   source_term: 'dada_maxEE1 | dada_maxEE2'
chimera_check_method: 'DADA2, removeBimeraDenovo, consensus'
chimera_check_param: not applicable
otu_clust_tool: 'DADA2, pool="pseudo"'
otu_clust_cutoff: 100
min_reads_cutoff: 
   default: 2
   options: 'unfiltered = 1 | filtered-trusted = 1 | filtered-analytic >= 2'
   source_file: decontam_workflow_config
   source_term: TBD (n_ton_removal)
min_reads_cutoff_unit: reads
min_reads_tool: 'DADA2 | decontam_workflow'
otu_db:
   - REVAMP:
      default: NCBI GenBank nt database, downloaded {nt_database_version}'
      source_file: REVAMP_config
      source_term: nt_database_version
   - SILVAngs:
      default: 'non-redundant SILVA SSU ref dataset, release {silva_db_release}'
      source_file: REVAMP_config
      source_term: silva_db_release
   - scikit-learn-silva:
      default: 'TBD'
      source_file: TBD
      source_term: TBD
otu_db_custom: not applicable
tax_assign_cat:
   - REVAMP:
      default: sequence similarity
   - SILVAngs:
      default: sequence similarity
   - scikit-learn-silva:
      default: probabilistic
otu_seq_comp_appr:
   - REVAMP:
      default: blastn >2.14.1+
   - SILVAngs:
      default: blastn 2.11.0+
   - scikit-learn-silva:
      default: 'no alignment, k-mer based'
tax_class_id_cutoff:
   - REVAMP:
      default: 60
      options: 'species = 97 | genus = 95 | family = 90 | order = 80 | class = 70 | phylum = 60'
      source_file: REVAMP_config
      source_term: taxonomyConfidenceCutoffs
   - SILVAngs:
      default: 86
   - scikit-learn-silva:
      default: not applicable
tax_class_query_cutoff:
   - REVAMP:
      default: 90
      source_file: REVAMP_config
      source_term: blastQueryCovCutoff
   - SILVAngs:
      default: 86
   - scikit-learn-silva:
      default: not applicable
tax_class_collapse:
   - REVAMP:
      default: 'Taxonomic levels were dropped to the lowest common ancestor (LCA), and were further dropped depending on % identity thresholds ({taxonomyConfidenceCutoffs}, see REVAMP readme).'
      source_file: REVAMP_config
      source_term: taxonomyConfidenceCutoffs
   - SILVAngs:
      default: 'No taxonomic levels are dropped from best BLAST hit matches due to high quality comprehensive nature of the reference dataset and taxonomy.'
   - scikit-learn-silva:
      default: TBD check - 'Starting at the highest level, the software checks % confidence against a set threshold. Once the confidence drops below that threshold, everything beyond that level is unassigned.'
tax_class_other: not applicable
screen_geograph_method: not applicable
screen_contam_0_1:
   - unfiltered:
      default: 0
   - filtered-trusted:
      default: 1
   - filtered-analytic:
      default: 1
screen_contam_method: 'TBD check - 1) The composition of the positive control is used to estimate a maximum vector contamination, which is then subtracted proportionally from all ASVs in the run to remove background tag jumping. 2) Next negative control contaminants are removed either as a wholesale removal of the impacted ASV or as a partial removal. 3) ASVs assigned to common contaminants are removed: human, food products, pets, common lab/consumable contaminants, laboratory controls.'
screen_nontarget_method: 'TBD check - 1) Contaminating off-target organisms were removed as described in screen_contam_method, removing common contaminants including human, food products, pets, common lab/consumable contaminants, laboratory controls.'
screen_other:
   - filtered-trusted:
      default: 'TBD check - In addition to the screening in screen_contam_method and screen_nontarget_method: 1) Remove positive and negative control samples. 2) Extreme low read depth sample removal. 3) Extreme low diversity sample removal (>99.9% of reads in one ASV). 4) Sample removal due to replicate dissimilarity distance from centroid above threshold.'
   - filtered-analytic:
      default: 'TBD check - In addition to the screening in screen_contam_method and screen_nontarget_method: 1) Remove positive and negative control samples. 2) Extreme low read depth sample removal. 3) Extreme low diversity sample removal (>99.9% of reads in one ASV). 4) Sample removal due to replicate dissimilarity distance from centroid above threshold. 5) Removal of singleton ASVs. 6) Removal of ASVs found in only one sample (no pattern of presence). Optional: 7) Remove n-ton ASVs (ASVs w/ less than n reads). 8) Low diversity sample removal. 9) Removal of unknown ASVs. 10) Low read depth sample removal. Exact application indicated in manuscript and manuscript code repository.'
otu_raw_description: 'No filtering outside of DADA2 default ASV denoising'
otu_final_description: this_DOI (link to decontamination screening section)
bioinfo_method_additional: this_DOI (paste link when published)

# NOAA PMEL Ocean Molecular Ecology terms
asv_method: 
   default: dada2pe
   options: dada2pe | dada2se
   source_file: REVAMP_config
   source_term: asv_method
discard_untrimmed: 1
dada2_trunc_len_f: not applicable
dada2pe_trunc_len_r: not applicable
dada2_trim_left_f: 
   default: 0
   source_file: REVAMP_config
   source_term: dada_trimLeft
dada2pe_trim_left_r:
   source_file: REVAMP_config
   source_term: dada_trimLeft
dada2_trim_right_f: 
   source_file: REVAMP_config
   source_term: dada_trimRight
dada2pe_trim_right_r:
   source_file: REVAMP_config
   source_term: dada_trimRight
dada2_max_ee_f:
   source_file: REVAMP_config
   source_term: dada_maxEE1
dada2pe_max_ee_r:
   source_file: REVAMP_config
   source_term: dada_maxEE2
dada2_trunc_q:
   source_file: REVAMP_config
   source_term: dada_trunQ
dada2_pooling_method: pseudo
dada2_chimera_method: 'removeBimeraDenovo, consensus'
dada2_min_fold_parent_over_abundance: not applicable
dada2_n_reads_learn:
   default: 1e11
   calculation: '{systemmemoryMB}*0.7*450000'
   source_file: REVAMP_config
   source_term: systemmemoryMB
   
---

# NOAA PMEL OME Bioinformatics Metabarcoding Protocol

## PROTOCOL INFORMATION

### Minimum Information about an Omics Protocol (MIOP)

- MIOP terms are listed in the YAML frontmatter of this page.
- See <https://github.com/BeBOP-OBON/miop/blob/main/model/schema/terms.yaml> for the list and definitions.

### Authors

- All authors known to have contributed to the preparation of this protocol should be listed, including those who filled in the template.
- Visit <https://orcid.org/> to register for an ORCID.

| PREPARED BY | AFFILIATION | ORCID | DATE |
| ------------- | ------------- | ------------- | ------------- |
| Sean McAllister | Ocean Molecular Ecology, NOAA PMEL & UW CICOES | <https://orcid.org/0000-0001-6654-3495> | 2025-05-08 |
| Samantha Setta | Ocean Molecular Ecology, NOAA PMEL & UW CICOES | <https://orcid.org/0000-0001-9075-7573> | 2025-04-30 |
| Shannon Brown | Ocean Molecular Ecology, NOAA PMEL & UW CICOES | <https://orcid.org/0000-0001-9808-2638> | 2025-04-30 |
| Zachary Gold | Ocean Molecular Ecology, NOAA PMEL | <https://orcid.org/0000-0003-0490-7630> | 2025-04-30 |

### Protocol Revision Record

- Version numbers start at 1.0.0 when the protocol is first completed and will increase when changes that impact the outcome of the procedure are made (patches: 1.0.1; minor changes: 1.1.0; major changes: 2.0.0).
- Release date is the date when a given protocol version was finalised.
- Description of revisions includes a brief description of what was changed relative to the previous version.

| VERSION | RELEASE DATE | DESCRIPTION OF REVISIONS |
| ------------- | ------------- | ------------- |
| 1.0.0 | 2025-MM-DD | Initial release |

### Related Protocols

- This section contains protocols that should be known to users of this protocol.
- Include the link to each protocol.
- Include the version number and release date (if available).
- Internal/External: "Internal" are derivative or altered protocols, or other protocols in this workflow. "External" are protocols from manufacturers or other groups.

| PROTOCOL NAME | LINK         | VERSION      | RELEASE DATE | INTERNAL/EXTERNAL |
| ------------- | ------------ | ------------ | ------------ | ----------------- |
| [decontam](https://github.com/benjjneb/decontam) | [Davis, Proctor, Holmes, Relman, Callahan 2018](https://doi.org/10.1186/s40168-018-0605-2)  | 2025-12-19  | External |
| [microDecon](https://github.com/donaldtmcknight/microDecon)  | McKnight 2021  | 2025-12-19  | External |
| [Gruinard decon](https://github.com/zjgold/gruinard_decon) | Gold 2020 | 2025-12-19 | External |

### Acronyms and Abbreviations

| ACRONYM / ABBREVIATION | DEFINITION |
| ------------- | ------------- |
| 12S | 12S ribosomal nucleic acid sequencing gene region |
| 16Sv4 | 16S ribosomal nucleic acid V4 gene region |
| 18Sv4 | 18S ribosomal nucleic acid V4 gene region |
| 18Sv9 | 18S ribosomal nucleic acid V9 gene region |
| ASV | Amplicon Sequencing Variant |
| BLAST | Basic Local Alignment Search Tool |
| CICOES | Cooperative Institute for Climate, Ocean, and Ecosystem Studies |
| COI | Cytochrome c oxidase subunit I gene region 
| DwC| Darwin Core |
| GBIF | Global Biodiversity Information Facility |
| ITS1 | Internal Transcribed Spacer 1 region |
| IUPAC | Content Cell |
| NCBI | National Center for Biotechnology Information |
| NOAA | National Oceanographic and Atmospheric Administration |
| OBIS | Ocean Biodiversity Information System | 
| ODE | Ocean DNA Explorer | 
|OME| Ocean Molecular Ecology|
| PMEL | Pacific Marine Environmental Laboratory |
| REVAMP | Rapid Exploration and Visualization through an Automated Metabarcoding Pipeline |
| SRA | 
| UW | University of Washington |

### GLOSSARY

| SPECIALISED TERM | DEFINITION |
| ------------- | ------------- |
| Denoising | The process of identifying and correcting sequencing errors to resolve true biological sequences (ASVs) from noise. |
| Demultiplexing| The process of sorting raw sequencing reads into sample-specific files based on unique index sequences (barcodes). |
| Index hopping | A phenomenon where indices are incorrectly assigned to the wrong library during sequencing, often mitigated by the UDIs (Unique Dual Indexes) mentioned in your protocol. |
| Negative control | Samples containing no intentional template DNA, used to monitor contamination during collection or processing.|
| Positive control  | A sample with a known DNA composition used to verify the success of the sequencing run and bioinformatic pipeline. |
| Singleton | An ASV that appears only once in the entire dataset (often filtered out to reduce noise). |
| SILVA  | Content Cell  |

## BACKGROUND

This document describes the required protocol to conduct insert name of the method/protocol.

### Summary

Insert a short description of the background for the method/protocol (e.g. why and for which purpose do you perform water sampling). Please provide a brief summary of your method including, as appropriate, a brief description of what techniques your best practice is about, which ocean environments or regions it targets, the primary sensors covered, what type of data/measurements/observing platform it covers, limits to its applicability.

### Method Description and Rationale

Insert a short description of the functioning principal of the methodology used in the protocol (i.e. how does the method work?). Please note that this is different from the step-by-step description of the protocol procedure.
Insert a short statement explaining why the specific methodology used in the protocol has been selected (e.g. it is highly reproducible, highly accurate, procedures are easy to execute etc….).

Data is run through this standard operating procedure at the sequencing run level, including all samples from the run, regardless of source project. The purpose of this is to maintain reproducibility and operate DADA2 as intended, allowing for proper error learning.

### Personnel Required

Insert the number of technicians, data managers, and scientists required for the good execution of the procedure

### Safety

No hazards or specialized safety equipment are associated with this protocol.

### Training Requirements

Specify the technical training required for the good execution of the procedure.

### Time Needed to Execute the Procedure

Specify how much time is necessary to execute the procedure.

## EQUIPMENT

| NAME | VERSION OR MODEL | MANUFACTURER OR CREATOR | REMARKS |
| ------------- | ------------- | ------------- | ------------- |
| Equipment |
| Laptop | Content Cell | Content Cell | e.g. needs at least 16 GB of RAM |
| Content Cell | Content Cell | Content Cell | Content Cell |
| Software |
| Content Cell | Content Cell | Content Cell | Content Cell |
| Content Cell | Content Cell | Content Cell | Content Cell |
| Code |
| Please include the links to the code you used for this analysis |
| e.g. link to the released version of a github repository  | Content Cell | Content Cell | Content Cell |
| Content Cell | Content Cell | Content Cell | Content Cell 

## STANDARD OPERATING PROCEDURE

### Raw Data Download and QA/QC

Raw reads (fastq.gz) are provided by the sequencing center, demultiplexed by sample and marker (see sequencing center BeBops listed in RELATED PROTOCOLS). Demultiplexing is done in two stages, depending on whether or not any sample•marker pair shares the same Illumina index. If all sample•marker pairs have a unique index, then demultiplexing is completed using pheniqs v2.1.0 only, with 0 allowed mismatches. If samples have unique indices, but markers do not, then first samples are separated with pheniqs v2.1.0 (0 mismatch), followed by marker separation using Cutadapt v3.4 (loose 30% allowed mismatch to primer; primer kept for a second more strict matching later in the workflow). On downloading reads to the local compute infrastructure, the fidelity of the downloaded file is checked via a `md5sum` check. Md5s are either supplied by the sequencing center with the original download (fastq.gz.md5) or are provided separately.

To check the quality of the sequencing run, [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is run on all fastq files (`~/path/to/FastQC/fastqc -o fastqc -t 10 *gz`). In addition for checking the run for a number of quality statistics, the "Per base sequence quality" is used to determine the best cutoff for the `dada_trimRight` (number of bases to trim from right or end of sequence) and `dada_trimLeft` (number of bases to trim from the left or start of sequence) parameters for the REVAMP configuration file.

### Rapid Exploration and Visualization through an Automated Metabarcoding Pipeline (REVAMP)

[REVAMP](https://github.com/McAllister-NOAA/REVAMP) ([McAllister et al., 2023](https://doi.org/10.5670/oceanog.2023.231)) is a published metabarcoding pipeline that integrates several tools for read processing, amplicon sequence variant (ASV) assignment, taxonomic assignment, and basic data visualization. The front end of REVAMP, including read trimming and ASV assignment, is used for all markers to produce a consistent set of ASVs to work from for downstream analysis. Taxonomic assignment in REVAMP is also used for downstream analysis on all markers, though some have additional alternative taxonomic classification methods (detailed below).

#### Stage and Prepare Data

Usually, data comes from the sequencing facility with facility sample tags and sequencing lane tag information that is cumbersome for downstream processing. OME sorts our data for each sequencing run into marker-specific folders, and renames them to OME eDNA sample numbers (E numbers) for use in the REVAMP workflow and downstream metadata merging and analysis.

At a bare minimum, REVAMP requires as input:
 * [sample metadata file](https://github.com/McAllister-NOAA/REVAMP?tab=readme-ov-file#sample-metadata-file--s) with a `Sample` header, though it is useful for downstream figure generation to provide other information from controlled vocabulary (i.e., `controls`, `sites`, `replicates`, `lat`, `long`, `group1`, `group2`) or free text "chemistry" entries.
 * [pipeline configuration file](https://github.com/McAllister-NOAA/REVAMP?tab=readme-ov-file#pipeline-configuration-file--p) with the following controlled vocabulary:
      * `primerF`: IUPAC nucleotide code for forward primer
      * `primerR`: IUPAC nucleotide code for reverse primer
      * `blastQueryCovCutoff`: Query coverage cutoff for BLASTn. Recommended: 90.
      * `systemmemoryMB`: System memory limit in megabytes.
      * `locationNTdatabase`: Path to the folder containing the prepared nt database. See next section.
      * `taxonomyCutoffs`: Percent identity confidence cut offs for assigning taxonomy to different taxonomic levels, from species to phylum.
         * Recommended for ribosomal RNA genes: 97,95,90,80,70,60
         * Recommended for protein encoding genes: 95,92,87,77,67,60
      * `dada_minlength`: Minimum post-trimming length passing DADA2.
      * `dada_phix`: Whether an internal PhiX sequencing control is present in the data (TRUE/FALSE).
      * `dada_trunQ`: Setting for quality score trimming of reads.
      * `dada_maxEE1`: Maximum number of expected errors in forward read.
      * `dada_maxEE2`: Maximum number of expected errors in reverse read.
      * `dada_trimRight`: Number of bp to trim from the right of the reads. Highly variable and based on sequencing run quality. Check FastQC results.
      * `dada_trimLeft`: Number of bp to trim from the left of the reads. Usually unnecessary (i.e., 0), but check FastQC results.
      * `blastMode`: Either `allIN`, `mostEnvOUT`, or `allEnvOUT`, referring to entries in the nt database that are kept or discarded when labelled with controlled unknown/unclassified vocabulary (see [ncbi_db_cleanup.sh](https://github.com/McAllister-NOAA/REVAMP/blob/main/ncbi_db_cleanup.sh)). Recommended: `mostEnvOUT`.
 * [figure configuration file](https://github.com/McAllister-NOAA/REVAMP?tab=readme-ov-file#figure-configuration-file--f):
    * See REVAMP readme. Not pertinent to bioinformatic processing.

REVAMP uses BLASTn against NCBI's nt database for taxonomic assignment. To run BLASTn, the nt database must be downloaded and prepared for use with REVAMP. This can be done using the shell script [ncbi_db_cleanup.sh](https://github.com/McAllister-NOAA/REVAMP/blob/main/ncbi_db_cleanup.sh) or by running each individual step from the script in the command terminal. Since the nt database is not version controlled, it is imperitive that the date of download is recorded to know the compatibility and relationship of different BLASTn runs. OME uses the same nt database download for the BLASTn for all markers on at least the same sequencing run.

#### Running REVAMP

Instructions for running REVAMP are included in the REVAMP [readme](https://github.com/McAllister-NOAA/REVAMP/blob/main/README.md).

OME runs REVAMP on all sequencing runs using the following command:

```
revamp.sh -p config.txt -f figure_config_file.txt -s sample_metadata.txt -r reads_dirs/marker -o marker_REVAMP -t 40 -y -g
```
Note that `-y` skips all intermediate prompts (choosing to automatically fill gaps in the taxonkit output between known levels). `-g` skips figure generation.

The purpose of REVAMP at the front of this workflow is to provide both a common set of ASVs for downstream analysis and to have the ability to apply at least one taxonomic method across all markers. Files used for downstream analysis include: ASV fasta file, ASV count table, ASV taxonomy assignment.

#### REVAMP Troubleshooting

When a sample is completely filtered out (zero remaining reads) at the Cutadapt, DADA2 part 1 (trim and filter), or DADA2 part 2 (learning error, dereplication, merge, ASV generation) steps (see the `run.log` for stats on each step), then REVAMP will fail ungracefully on subsequent steps. The best solution is to remove the offending sample from the run entirely (remove from raw_read folder and sample_metadata file) and start the pipeline again from scratch.

### Taxonomic Classification

#### `REVAMP`

Citation for the Rapid Exploration and Visualization through an Automated Metabarcoding Pipeline (REVAMP) tool: [https://doi.org/10.5670/oceanog.2023.231](https://doi.org/10.5670/oceanog.2023.231).

Modified from the REVAMP [workflow](https://github.com/McAllister-NOAA/REVAMP/tree/main?tab=readme-ov-file#revamp-workflow): 
* ASVs are blasted against the NCBI `nt` database (BLASTn; `subject_besthit`; `max_target_seqs 4000`) ([Camacho et al., 2009](https://doi.org/10.1186/1471-2105-10-421)), exporting percent identity, length of hit, subject taxonomy IDs (taxIDs), and subject accession in tab-delimited format.
* After the tab-delimited BLAST output file is created, the file is then reformatted to simplify the results to only include taxIDs from all best percent identity matches longer than a user-supplied query coverage cutoff (90% by length, in this case). Taxonomy assignment with REVAMP is intentionally conservative given: 1) the reference databases are incomplete for some taxa and heavily sampled for others, and 2) markers vary in their ability to resolve different taxa to species (e.g., [Gold et al., 2021](https://doi.org/10.1111/1755-0998.13450)).
* Before assessing taxonomies, TaxonKit (v 0.5.0) ([Shen and Ren, 2021](https://doi.org/10.1016/j.jgg.2021.03.006)) is run on all discovered taxIDs in the reformatted BLASTn file, with the results reformatted to include kingdom, phylum, class, order, family, genus, and species (`K/P/C/O/F/G/S`) assignments only.
* In cases where an intermediate taxon is missing in this string, REVAMP will automatically fill from a lower taxonomy assignment: for example, a `K/P/-/O/F/G/S` string would be replaced with `K/P/O__c/O/F/G/S` to indicate that the class of interest contains the order below it.
* In REVAMP, taxonomy is then assigned by merging the taxonomic hierarchy of all best BLASTn taxID hits to the lowest common ancestor (i.e., deepest identical taxonomic assignment). This means, for example, that if the best BLASTn hits match several species from the same genus, the ASV will only be identified to the genus level.
* The only other factor influencing the assignment is the confidence of the taxonomic depth as determined by percent identity cutoffs for each taxonomic level. While different clades vary widely in this respect ([Gebhardt and Knebelsberger, 2015](https://doi.org/10.1007/s10152-015-0434-7); [Zhang and Bu, 2022](https://doi.org/10.3390/insects13050425)), users can choose a set of taxonomic cutoffs per marker gene. For this marker, we use the recommended confidence cut offs for rRNA genes (97,95,90,80,70,60), which means that an ASV with best blast hit percent identity ≥97% can be assigned to the species level, 97%>x≥95% to genus, 95%>x≥90% to family, 90%>x≥80% to order, 80%>x≥70% to class, 70%>x≥60% to phylum, and <60% is assigned to unknown. Note that this confidence-based taxonomic trimming is completed after the last common ancestor analysis trimming.
* Taxonomy is output as a tab-delimited file with one ASV per line, followed by the seven levels of taxonomy assignment. Unassigned levels are designated as an `NA`.
* Summary:
   * BLAST-based best percent identity hit (over 90% query coverage by length)
   * Lowest common ancestor of all best hits
   * Confidence thresholds set the possible depth of assignment based on the percent identity quality of the hit
   * Taxonomic assignments are based on the NCBI taxonomic hierarchy, simplified to seven levels only (`K/P/C/O/F/G/S`)
   * Results will change depending on the contents of the NCBI `nt` database, which is why the date of download should be noted

#### `SILVAngs`

Citation for the SILVAngs tool: [https://doi.org/10.1016/j.jbiotec.2017.06.1198](https://doi.org/10.1016/j.jbiotec.2017.06.1198). The tool is a web-based service found [here](https://ngs.arb-silva.de/silvangs/). Detailed description of the data processing pipeline is described in the [User Manual](https://www.arb-silva.de/fileadmin/silva_databases/sngs/SILVAngs_User_Guide.pdf). The [arb-silva](https://www.arb-silva.de/) team is well-known and trusted in the realm of microbial taxonomic classification and reference database management, having maintained a curated database of the small and large subunit ribosomal RNA genes for more than 30 years. The SILVAngs tool allows a user to submit sequences to be automatically classified against the curated SILVA taxonomy backbone to the genus level.

A summary of the SILVAngs pipeline and OME's use of it:
* The `ASVs.fa` amplicon sequence variant output from the REVAMP pipeline is submitted to the SILVAngs platform on a per sequencing run, per marker basis.
* Default analysis settings are used with the exceptions: `sequence identity` set to `1` (these already represent ASVs, aka 100% operational taxonomic units).
   * SILVAngs, as run on default settings, runs BLASTn against the SILVA database (release 138.2), and only assigns reads where the average of the percent identity and alignment coverage are greater than 93%. No common ancestor trimming is done due to the nature of the non-redundant curated database.
* The project is executed by the SILVAngs pipeline, and the result archive is downloaded.
* To create a similar data product to REVAMP, which includes ASV counts in a biom format and ASV taxonomy in a separate file, we then use REVAMP to reformat the SILVAngs output.
   * Note: Despite setting `sequence identity` set to `1`, there are still rare instances where multiple ASVs can be clustered together by the SILVAngs pipeline. REVAMP will look at the clustering file to copy the taxonomic assignment of the cluster reference to other ASVs in that cluster.
* The original REVAMP output folder that was used to create the original ASV assignments is copied to a new location. Because of the checkpoint system in REVAMP, we can overwrite the original REVAMP taxonomy of this copied folder with the new SILVAngs taxonomy by:
   * Deleting all the checkpoint lines in `progress.txt` after `dada2_Finished=TRUE`.
   * Rerunning REVAMP with the `-e` flag and pointing to the config files internal to the output folder.
   * E.G. `revamp.sh -p Parada16S_REVAMP_out/config_file.txt -f Parada16S_REVAMP_out/figure_config_file.txt -s Parada16S_REVAMP_out/sample_metadata.txt -r . -o Parada16S_REVAMP_out -t 6 -y -g -e`
* REVAMP will skip the BLASTn steps and instead prompt the user for the SILVAngs output and SILVA taxonomy files:
   * `Enter the location of the SILVAngs ssu or lsu results directory (i.e. ~/Downloads/results/ssu)`. This is simply the path to the results archive `ssu` folder.
   * `Enter the location of the reference taxonomy map for current SILVA database: i.e., tax_slv_ssu_138.1.txt`. This file can be downloaded from arb-silva [here](https://www.arb-silva.de/current-release/Exports/taxonomy).
* REVAMP will output the necessary ASV taxonomy file, based on the SILVA taxonomy hierarchy, and simplifyit  to seven levels only (`K/P/C/O/F/G/S`). Note: All species-level assignments in this file are `NA`.

##### Troubleshooting
There is a rare designation in the SILVAngs output that sets the taxonomy assignment to `silva||0|`. This entry does not exist in the SILVA taxonomy database, and is meant to indicate an unassigned/unknown hit. However, REVAMP cannot yet deal with this issue (user will see perl error `Use of uninitialized value`), and it is necessary to replace `silva||0|` with `Unknown` in the `results/ssu/exports/x---ssu---otus.csv` file. Note that the assignment `ncbi||0|` in the same line should also be replaced with `Unknown`.

#### `scikit-learn-silva`

The `ASVs.fa` amplicon sequence variant output from the REVAMP pipeline 16Sv45 rRNA region is classified using [Qiime2 (v.2024.10)](https://qiime2.org/) (Bolyen et al., 2019) feature classifier's [naive bayesian classifier scikit-learn](https://scikit-learn.org/stable/modules/naive_bayes.html). Sci-kit learn classifiers is trained using the [silva (v138.1)](https://docs.qiime2.org/2022.11/data-resources/) for 16Sv4 rRNA (Quast et al., 2013). Each database is curated, extracted to the region of interest, and used to train a naive Bayesian taxonomic classifier by:

1) Downloading the full-length silva (v138.1) taxonomy and sequences.
```
wget https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza
```

2) Extracting the region of interest from all sequences in the database with primers used for amplification (see [NOAA PMEL OME Github](https://github.com/NOAA-PMEL/Ocean-Molecular-Ecology)).

```
# extract region of interest with primers and target length
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer CCGYCAATTYMTTTRAGTTT \
  --p-min-length 300 \
  --p-max-length 500 \
  --o-reads silva_v138.99_v4_extracts.qza
# --p-min-length and --p-max-length vary by region and are included in the config file.

# Get taxa that remain from extracts
qiime rescript filter-taxa \
    --i-taxonomy silva-138-99-tax.qza \
    --m-ids-to-keep-file silva_v138.99_v4_extracts.qza \
    --o-filtered-taxonomy silva_v138.99_v4_tax.qza

# number of seqs, lengths, etc.
qiime feature-table tabulate-seqs \
  --i-data silva_v138.99_v4_extracts.qza \
  --o-visualization silva_v138.99_v4_extracts.qzv

```
3) Adding sequences that contain the region of interest but do not contain the primer region by alignment to the extracted sequence database. This step is done twice for some metabarcoding regions until ~50-70% of query sequences are retained following [recommendations](https://forum.qiime2.org/t/using-rescripts-extract-seq-segments-to-extract-reference-sequences-without-pcr-primer-pairs/23618) for qiime2 scikit-learn classifiers (see config file).

```
# add sequences that did not contain primers, but do have the region of interest
qiime rescript extract-seq-segments \
    --i-input-sequences silva-138-99-seqs.qza \
    --i-reference-segment-sequences silva_v138.99_v4_extracts.qza \
    --p-perc-identity 0.8 \
    --p-min-seq-len 150 \
    --p-threads 5 \
    --o-extracted-sequence-segments silva_v138.99_v4_extracts-01.qza \
    --o-unmatched-sequences silva_v138.99_v4_unmatched-extracts-01.qza \
    --verbose 
# --p-perc-identity and --p-min-seq-len vary by region and are included in the config file.

# Compare counts of seqs
qiime feature-table tabulate-seqs \
  --i-data /home/poseidon/setta/ref_db/silva/silva_v138.1/silva_v138.99_v4_extracts-01.qza \
  --o-visualization /home/poseidon/setta/ref_db/silva/silva_v138.1/silva_v138.99_v4_extracts-01.qzv
```

*Note: Step 3 is repeated until ~50-70% query sequnces retained, with the `--o-extracted-sequence-segments` output file used in the `--i-reference-segment-sequences` input for the next iteration. After the first iteration, increase the `--p-perc-identity` to 0.9. For the 16Sv4 region, only needed 1 iteration.*

4) Filtering out taxonomies for sequences removed during the extraction of full length database.
```
# Get taxa that remain:
qiime rescript filter-taxa \
    --i-taxonomy silva-138-99-tax.qza \
    --m-ids-to-keep-file silva_v138.99_v4_extracts-01.qza \
    --o-filtered-taxonomy silva_v138.99_v4_tax-01.qza
# replace `--m-ids-to-keep-file` with extracted sequence file from the last iteration of step 3.
```

5) Dereplicating the extracted reads, retaining identical sequence records that have differing taxonomies. 

```
qiime rescript dereplicate \
    --i-sequences silva_v138.99_v4_extracts-01.qza \
    --i-taxa silva_v138.99_v4_tax-01.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva_v138.99_v4_uniq-seqs.qza \
    --o-dereplicated-taxa silva_v138.99_v4_uniq-tax.qza
```

6) Evaluating reference taxonomy and sequences.

```
# original databases
qiime metadata tabulate \
    --m-input-file silva-138-99-seqs.qza \
    --o-visualization silva-138-99-seqs.qzv
# extracted and dereplicated databases  
qiime metadata tabulate \
    --m-input-file silva_v138.99_v4_uniq-seqs.qza \
    --o-visualization silva_v138.99_v4_uniq-seqs.qzv
# extracted and dereplicated databases    
qiime rescript evaluate-taxonomy \
    --i-taxonomies silva_v138.99_v4_uniq-tax.qza \
    --o-taxonomy-stats silva_v138.99_v4_uniq-tax.qzv
```

7) Training the taxonomic classifier for each region of interest using qiime2's feature-classifier.
```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva_v138.99_v4_uniq-seqs.qza \
  --i-reference-taxonomy silva_v138.99_v4_uniq-tax.qza \
  --o-classifier silva_v138.99_v4_uniq-classifier.qza
```

8) The trained classifiers are used to classify ASVs using qiime2's feature classifier.
```
# import fasta file from revamp into qiime formatted file
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 16Sv4_ASVs.fa  \
  --output-path 16Sv4_ASVs.qza

# run classifier
qiime feature-classifier classify-sklearn \
   --i-reads 16Sv4_ASVs.qza \
   --i-classifier silva_v138.99_v4_uniq-classifier.qza \
   --p-n-jobs 4 \
   --o-classification 16Sv4_tax_silva_sklearn.qza

# export taxonomy into tab-delimited file
qiime tools export \
  --input-path 16Sv4_tax_silva_sklearn.qza \
  --output-path 16Sv4_tax_silva_sklearn
```

The final output file is the taxonomy.tsv file within the 16Sv4_tax_silva_sklearn directory with `Feature ID` (ASVs), `Taxon` (taxonomy separated by semicolons), and `Confidence` (confidence score of the lowest taxonomic level). Silva taxonomy output includes 7 levels: domain, phylum, class, order, family, genus, species.

#### Taxonomy Choices

After running all of the above taxonomic classifiers, we still need a choice taxonomy to submit to Ocean DNA Explorer ([ODE](https://www.oceandnaexplorer.org/)) and Ocean Biodiversity Information System ([OBIS](https://obis.org/)). This choice differs between markers and is currently under investigation (Pers. comm. Sam Setta, Zack Gold, Sean McAllister). Once that choice is made, it is then crosswalked using [edna2obis](https://github.com/aomlomics/edna2obis) with the World Register of Marine Species ([WoRMS](https://www.marinespecies.org/)) and the OBIS and Global Biodiversity Information Facility ([GBIF](https://www.gbif.org/)) taxonomic backbones for submission. This field represents the best taxonomic assignment (given for the controlled Darwin core observation records) for each ASV, though independent assignments through the methods listed above are still kept as optional user fields. 

### Data Decontamination and Quality Assurance Processing Workflow
![Sequencing-data-decontamination-and-quality-assurance (1)](https://github.com/user-attachments/assets/7a46aa0d-1247-4550-b19a-0989b80c4af4)

Unfiltered ASV tables output from any bioinformatics pipeline have the potential for a number of artefacts and issues that should be dealt with before providing those results to the public. OME has chosen to provide ASV tables at three stages of analysis: 
   1) unfiltered - The original ASV table with no additional filtering. Submitted to Ocean DNA Explorer ([ODE](https://www.oceandnaexplorer.org/)).
   2) filtered-trusted - An ASV table filtered to remove ASVs and samples with clear issues affecting their use in any context. Submitted to [ODE](https://www.oceandnaexplorer.org/) and Ocean Biodiversity Information System ([OBIS](https://obis.org/)).
   3) filtered-analytic - Starting from the filtered-trusted ASV table, additional filtering is applied for a given, more detailed analytic context. Submitted to [ODE](https://www.oceandnaexplorer.org/) and as part of a code repository on GitHub associated with a particular manuscript.

The OME workflow for data decontamination and quality assurance is available [here](https://github.com/DiatomSetta/Sequencing-data-decontamination-and-quality-assurance/tree/main). In brackets below are the types of data being filtered, either ASVs or Samples.

**Part 1 – Sequencing Run-Based Filtering**

**Step 1 – Filter Read Type** 

By default, REVAMP only produces ASV tables from merged forward and reverse reads. However, some pipelines produce multiple datasets from single-end and merged reads (e.g., [Anacapa](https://github.com/limey-bean/Anacapa)), and this step allows the user to filter those read types depending on use. Merged reads are more accurate, and OME primarily focuses on that result. However, there are applications of single-end ASVs, including 1) instances where amplicons are larger and provide no overlap of the forward and reverse reads. This is the case for Eukaryota with the `ssu16sv4v5_parada` assay, for instance. 2) Failed sequencing runs where only one read (forward or reverse) was successful.

**Step 2 – Filter Tag-jumping [ASV]**

If the sequencing run includes a high-quality positive control, that control can be used to estimate the occurrence of tag-jumping or index hopping. The ASV(s) associated with the positive control sample are first identified, after which the maximum relative proportion of environmental reads jumping into the positive control is calculated as maximum vector contamination. Next, that maximum contamination proportion is subtracted from each ASV's total read count.

**Step 3 – Filter Negative Control Contaminants [ASV]**

If a sequencing run has negative controls sequenced, those negative controls can be used to remove contaminants. It is important to consider the source of the negative control so that the filtering can be applied appropriately to each sample. For example, OME collects and sequences negative controls throughout cruises (field blanks) and during extraction (extraction blanks). Depending on which samples each of these negative controls is associated with, the clearance of contaminating sequences may not apply to all samples on a sequencing run evenly. For those samples that negative controls do apply to, OME filters using two strategies: 1) wholesale removal of the impacted ASV, or 2) partial removal.

For the partial removal, ASVs found in the negative control are slated for removal if the following conditions are met:
   1) ASV occurs in less than 10 samples.
   2) ASV does not occur as an abundant (≥10% per sample max occurrence) ASV in any sample.
   3) ASV abundance is outside a 95% confidence interval for a normal distribution from the mean total reads and less than the total sum of reads in the negative control.

Both approaches are valid. OME primarily uses the partial removal approach, since there is a high likelihood that an abundant ASV from real environmental samples might tag jump into the negative control. If it isn't fully removed during Step 2, then this would lead to the removal of an ecologically important ASV from the dataframe.

TBD - check wording.

**Step 4 – Remove Positive and Negative Control Samples [Sample]**

After they are used in steps 2 and 3, the positive and negative control samples can be removed from the dataset.

**Step 5 – Remove Known Contaminants and Off-targets [ASV]**

There are several typical contaminants that may make their way into environmental samples, particularly for low DNA applications of eDNA metabarcoding. These can include DNA from the humans processing the samples, food items (particularly meats), pets (dogs, cats, etc.), and lab controls (gblocks for positive controls, positive controls from other experiments). In addition, some markers were "designed" to target particular taxa, but may in fact amplify other organisms that were not originally intended as the target. This step allows users to designate known contaminants and off-targets to remove them from the data.

Users can supply a list of typical lab contaminant taxonomies (one per line), which should include the potential for incomplete taxonomic assignment of the contaminant in the current dataset. Users can also supply the sequence of contaminants for direct removal of matching ASVs (independent of taxonomy). This latter is useful for removing the sequences of synthetic controls, which may or may not be assigned to the appropriate taxonomy in a traditional workflow.

OME known contaminant list:
* Eukaryota;Chordata;Mammalia;Primates;Hominidae;Homo;Homo sapiens #Human contamination potential during collection and lab processing
* Eukaryota;Chordata;Mammalia;Primates;Hominidae;Homo;NA #Human contamination potential during collection and lab processing
* Eukaryota;Chordata;Mammalia;Primates;Hominidae;NA;NA #Human contamination potential during collection and lab processing
* Eukaryota;Chordata;Mammalia;Primates;NA;NA;NA #Human contamination potential during collection and lab processing, primates unlikely in the NE Pacific Ocean
* Eukaryota;Chordata;Mammalia;Artiodactyla;Bovidae;Bos;Bos taurus #Domestic cattle contamination potential from foods
* Eukaryota;Chordata;Mammalia;Artiodactyla;Suidae;Sus;Sus scrofa #Wild boar contamination potential from foods
* Eukaryota;Chordata;Aves;Galliformes;Phasianidae;Gallus;Gallus gallus #Chicken contamination potential from foods
* Eukaryota;Chordata;Mammalia;Carnivora;Canidae;Canis;Canis lupus #Wolf/dog contamination potential from pet dander/hair
* Eukaryota;Chordata;Mammalia;Carnivora;Felidae;Felis;Felis catus #Cat contamination potential from pet dander/hair
* Eukaryota;Chordata;Aves;Accipitriformes;Accipitridae;Harpagornis;Harpagornis moorei #Haast's eagle gblock positive control ssu12sv5v6_mifish_u_sales
* Eukaryota;Chordata;Mammalia;Artiodactyla;Camelidae;Camelus;Camelus dromedarius #Dromedary camel gblock positive control ssu16sv4v5_parada
* Eukaryota;Bigyra__p;Bigyra;Bicosoecida;Halocafeteria__f;Halocafeteria;Halocafeteria seosinensis #halophilic nanoflagellate gblock positive control ssu18sv9_amaralzettler & ssu18sv8_machida
* Eukaryota;Bacillariophyta;Bacillariophyceae;Naviculales;Diadesmidaceae;Luticola;Luticola ventricosa #Freshwater Antarctic diatom gblock positive control ssu18sv4_stoeck & ITS1_sterling
* Eukaryota;Chordata;Aves;Dinornithiformes;Dinornithidae;Dinornis;Dinornis giganteus #Giant Moa gblock positive control COI_1835-2198_lerayfolmer
* Eukaryota;Chordata;Mammalia;Artiodactyla;Delphinidae;Orcaella;Orcaella brevirostris #Irrawaddy river dolphin gblock positive control dLoop_Baker
* Eukaryota;Chordata;Aves;Columbiformes;Raphidae;Raphus;Raphus cucullatus #Dodo bird gblock positive control lsu16s_2434-2571_kelly
* Eukaryota;Chordata;Actinopteri;Perciformes;Sebastidae;Sebastes;Sebastes inermis #Japanese rockfish gblock positive control cytB_MiSebastes
* Eukaryota;Chordata;Mammalia;Perissodactyla;Rhinocerotidae;Coelodonta;Coelodonta antiquitatis #Woolly rhinoceros gblock positive control lsu16s_2051-2438_mideca
* Eukaryota;Chordata;Mammalia;Sirenia;Dugongidae;Hydrodamalis;Hydrodamalis gigas #Steller's sea cow gblock positive control COI_MollCOI253
* Eukaryota;Chordata;Mammalia;Pilosa;Megatheriidae;Megatherium;Megatherium americanum #Giant sloth gblock positive control lsu16s_Mol16S & DC_MA_nad2
* Eukaryota;Chordata;Aves;Columbiformes;Columbidae;Ectopistes;Ectopistes migratorius #Passenger pigeon gblock positive control Pycno_EM_nad5
* Eukaryota;Chordata;Mammalia;Artiodactyla;Camelidae;Camelus;NA #Camel gDNA positive control ssu16sv4v5_parada_OSUmod & ssu18sv8_machida_OSUmod
* Eukaryota;Chordata;Mammalia;Carnivora;Mustelidae;Mustela;Mustela putorius #Ferret gDNA positive control ssu16sv4v5_parada_OSUmod & ssu18sv8_machida_OSUmod
* ```GAGCGGAGAGATAACTATAATGCTGCTTATAATCTCTCTGTCGTTGCTAACGTTTATAGTCTAGTCTCATTATAATTGTATGCTATTGAGGCATTGATTAATACTGGAAAACATTTGAAATAAACTAATTTATACGACAGAAATCGTGTACCTACTAAATCTCTTTAATGTAAGTTCTGACTAATTCGTACTTTGTTAAGAACTTACATTTTAATAATAGAGGATATATGTTTTATTTTTATGATCTATTGATGTTCTTAAGACTGCAATTTATATAATGAGGTAATATTTGCGGTAAATCCTAGTGCAATGGCAATTTTTTACTTTTGTTCTAAAGAAGAGATAACGTGAGTGCAGTTATCATTAATGTAGAAATTGGAAAGATTCTTGGGCCTCCACCTTTAGATAGTGTTTACTCTTTTATAAAGGAGCTATTAATTATGTCTTGCGAAGATTCAAAAAGGTAAGTCAATTTGGCTGATTCGAAAAGACGGACTTCAAAGTTACCTAACAATAGTTGTGGGTCCGTAACAAAATCTTTTTAATAAGTCTCCGTAAGTGTTGGTTGAATAGCCCTGATCGGTTTAACCAGAAGTTCAAATATATTATTTTATGATTTTTGAGTCAGAATGTGTTACCTTACAGAAATTAAGATCG``` #Random 32% GC sequence gblock positive control COI_Folmer

Off-targets are specified by taxonomy or ASV length. For example, ssu12sv5v6_mifish_u_sales has a known bacterial off-target greater than 210 bp. For the MiFish assay, off-targets can be removed by size selection (e.g., keep <210 bp) or taxonomy (e.g., keep Eukaryota).

**Part 2 – Merging Datasets**

Data from a single project may end up on multiple sequencing runs. This part of the workflow allows us to combine multiple runs to make a single project-based dataset. ASVs with identical sequences are merged, new count tables generated, and conflicting taxonomiesare  identified and sorted.

**Part 3 – Depth, Diversity, Replicate Filtering**

**Step 6 – Remove Extreme Low Sequencing Depth Samples [Sample]**

Remove any sample with fewer than 1,000 total reads. These samples are likely poorly sequenced, and the data are not a good representation of the environment.

**Step 7 – Remove Extremely Low Diversity Samples [Sample]**

Remove any sample that contains a single ASV with >99% of the reads mapping to it. This situation is unlikely in the environment and likely indicates a problematic sample.

**Step 8 – Dissimilarity Between Replicates [Sample]**

For this analysis, sample relationships are graphed in ordination space, and a centroid is calculated for any pair of replicates. If a sample within a group of replicates has a distance to that centroid greater than a maximum distance cut off (TBD), then it is removed from the dataset.

**Filtered-Trusted Dataset**

At this point in the workflow, the filtered-trusted dataset is peeled off for submission to [ODE](https://www.oceandnaexplorer.org/) and [OBIS](https://obis.org/).

**Part 4 – Final Analytic Filtering And Optional Steps**

**Step 9 – Remove Singletons [ASV]**

Remove any ASV represented by only a single read across the dataset.

**Step 10 – Pattern of Presence [ASV]**

Remove any ASV that is only found in one sample in the dataset (or user-set sample count cutoff).

**Step 11 – Remove n-tons [ASV]**

The default for this step is to skip this filtering, though a user can set any total read count cutoff for removing ASVs from the dataset.

**Step 12 – Remove Low Diversity Samples [Sample]**

The default for this step is to skip this filtering, though a user can set any cutoff to eliminate samples that are majority dominated by a single species (x%).

**Step 13 – Remove ASVs Assigned to Unknowns [ASV]**

Remove any ASV assigned to a list of taxa strings that are associated with Unknown assignments. For REVAMP: `Unknown;NA;NA;NA;NA;NA;NA` or `Environmental Unknown;NA;NA;NA;NA;NA;NA`.

**Step 14 – Remove Low Sequencing Depth Samples [Sample]**

The default for this step is 10,000 reads. A user can supply any cutoff for removing samples based on their total read count/depth.

**Filtered-Analytic Dataset**

After filtering, which can be more or less stringent depending on needs, the analytic dataset is ready for submission to [ODE](https://www.oceandnaexplorer.org/) or for analysis for a manuscript. The final dataset and filtering choices will be supplied with the manuscript-specific GitHub code repository.
   
### Data and Metadata Submission Strategy

The OME datamanagement plan includes submission of all data to three distinct repositories:

1) Raw reads are submitted to NCBI's sequence read archive (SRA)
2) Unfiltered ASV tables and taxonomy assignments are submitted to [ODE](https://www.oceandnaexplorer.org/) and a manuscript-specific code repository.
3) Filtered-Trusted ASV tables and taxonomy assignments are submitted to [ODE](https://www.oceandnaexplorer.org/) and [OBIS](https://obis.org/).
4) Filtered-Analysic ASV tables and taxonomy assignments are submitted to [ODE](https://www.oceandnaexplorer.org/) and the manuscript-specific code repository.
5) Sample metadata are submitted with raw read records to NCBI's SRA, to ODE and OBIS through FAIR eDNA terms mapped to Darwin core controlled vocabulary, through the manuscript-specific code repository, and through an archive with NCEI.

### Marker-specific recommendations – ssu16sv4v5_parada / Universal-16S-V4V5-Parada

TBD

## GUIDE TO ARCHIVED METHODOLOGY

The contents of this archive should allow your analysis to be reproduced exactly as you intended it. This document provides guidance on the contents of each partner's compressed archive of in-silico methods. This document should be part of that same archive, serving as an extended README. Below, please find guidance on what this archive should include. When describing the contents of the archive, please give precise file names and relative paths to the files.

### Archive Content

To reproduce the in-silico analysis, please provide one of the following (in order of decreasing preference)

1. Jupyter, R notebook(s) or equivalents
2. Downloaded an archive of (the released version of) your GitHub repository
3. Individual scripts

In each of the above cases, guidance and documentation for all the steps you took to perform the in-silico analysis should be included. In case 1., code and documentation are integrated. In cases 2. and 3., in-line comments may be provided, however, these are not generally sufficient as documentation. In those cases, please provide a step-by-step protocol on how and when to run each script in the Execution Procedure section below.

Please include a script on **data acquisition** (e.g., documentation and code to pull sequences from INSDC, access sequences on an institutional FTP server, download metadata files, check file integrity via md5 checksum). Please add sufficient detail, so that the partners only have to install the software, run this script and will then have all the data needed to perform any analysis described below.

### Code
Here please describe each file containing code, including its purpose, its input, its output. Please provide the names and the relative paths to this documentation.

### Code documentation
Here, please indicate if your documentation is with the code (in a code notebook) or stored separately. In-line comments are not considered documentation. If the documentation is stored separately, please provide the names and the relative paths to this documentation.

### Metadata
Please provide link(s) to the files containing metadata about your sequence data (e.g. environmental data, procedural data). Please see the MIxS compliant metadata guidance.

Auxiliary files
e.g. mapping files, test/dummy files, colour palette

## Execution Procedure

Please fill out this section if you have not already documented it as part of your R, Jupyter, or similar notebook. In this section, please provide a step-by-step guidance on how and when to run each component of your code.

## Quality Control

In this section please include the names and paths that can be used to validate that operations were successful. If such checks were done during the execution procedures, please note this here. We recommend identifying such steps with in-line tags (e.g. “#QC”).

### Basic Troubleshooting Guide

Identify known issues associated with the procedure, if any.
Provide troubleshooting guidelines when available.

## REFERENCES

Bolyen, E., Rideout, J.R., Dillon, M.R. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852–857 (2019). [https://doi.org/10.1038/s41587-019-0209-9](https://doi.org/10.1038/s41587-019-0209-9)

Davis, N.M., Proctor, D.M., Holmes, S.P. Relman, D.A., & Callahan, B.J. (2018). Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome 6, 226. [https://doi.org/10.1186/s40168-018-0605-2](https://doi.org/10.1186/s40168-018-0605-2)

Guillou, L., Bachar, D., Audic, S., Bass, D., Berney, C., Bittner, L., ... & Christen, R. (2012). The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote small sub-unit rRNA sequences with curated taxonomy. Nucleic acids research, 41(D1), D597-D604. [https://doi.org/10.1093/nar/gks1160](https://doi.org/10.1093/nar/gks1160)

Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Opens external link in new windowNucl. Acids Res. 41 (D1): D590-D596.[https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)

DADA2 Callahan et al., 2016

Cutadapt (Martin, 2011)

REVAMP ([McAllister et al., 2023]

Link to any documents such as software guidelines, images, etc that support this protocol. Please include a short note describing the document's relevance.
