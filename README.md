# belmonte_lab

**Collection of scripts written for the Belmonte lab.**

### Formatting Programs:
- **`ExtractLines`**: Extract lines that contain a target string (old script, grep is much more efficient)
- **`PWM_to_IUPAC`**: Convert PWM (Position Weight Matrix) to IUPAc sequences (java)
- **`ProcessedGOFormatter`**: Process SeqEnrich GO term output into a more usable format, filters p-values (old script, could do much better in R)

---

### Methylation Project:
- **`2017_08_01_bsmap_methylKit.R`**: Analyse data from BSMAP methylation call files with methylKit (Also create BED4 DMR files)
- **`bedtools_cmds.txt`**: Bash commands to create gene flanks and intersect gene/flanks with DMRs (Differentially Methylated Regions)
- **`build_cytosine_bed4.R`**: Make BED4 files from our BSMAP methylation call files (after they went throught methratio.py)
- **`build_gene_expression_bed5.R`**: Make BED5 files of differential gene expression data. The information is taken from several files and put together in a useful format.
- **`correlate_meth_gene_expr.R`**: Older of the 2 correlation scripts, methylation level in a DMR vs gene expr the DMR overlaps with or flanks. Not very pretty or useful results-wise because it used older fpkm data and slightly messier input but could be useful as a reference.
- **`corr_mean_meth_vs_expr.R`**: Correlates the mean methylation ratio in a feature (gene or its flanks) with that feature's expression level (fpkms from cuffdiff). Makes scatterplots.
- **`extract_filter_intersectBed.R`**: Process the output of bedtools_cmds into a more usable format. Includes cmds for SeqEnrich prep
- **`methylation_heatmap_and_bargraph.R`**: Visualize methylation data by context, sample, and chromosome, visualize transposable element, gene, and CpG island density. (Mostly bar charts currently, heatmaps are not turning out great)

---

### Pipeline:
- **`SeqAnalysisPipeline`**: (perl) Guides user through the Belmonte lab's Sequencing Analysis Pipeline (Trimmomatic, TopHat2 with Bowtie2-build, Cuffquant, Cuffnorm, and CuffDiff) using cmd line prompts. Does preliminary error checking, saves a copy of the commands that were run, and organizes the output. Intructions/tutorial PDF file included in repository.


---

### siRNA Project:
- **`analysis\_of\_siRNA\_meth\_intersect.R`**: Compute the average methylation ratio per siRNA locus from Shortstack output based on the output of intersectBed using the bed files of both data sets prepared in this project.
- **`blastn\_print\_top\_hits.py`**: Do a local blastn of a given query file and print the top hit for each query sequence. Essentially a python wrapper and parser. Saves the full blast results in a .xml file as a byproduct.
- **`build_cytosine_bed4\_siRNA\_version.R`**: Edited version of the same script in the **Methylation Project**. Has a lower filtering threshold and updated input file locations.
- **`intersectBed\_cmds.txt`**: Documentation of commands used to intersect the methylation and ShortStack siRNA bed files.
- **`shortstack\_blastn\_to\_bed4.py`**: Create bed4 files for intersection from the .xml blastn outoput file of ShortStack siRNAs.
- **`shortstack\_cmds.txt`**: Documentation of commands used for ShortStack analysis.
- **`standalone_blastn_adventures.txt`**: Documentation of commands used for setting up and running blastn locally.
