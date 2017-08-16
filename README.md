# belmonte_lab

Collection of scripts written for the Belmonte lab.

Formatting Programs:
- ExtractLines: Extract lines that contain a target string (old script, grep is much more efficient)
- PWM_to_IUPAC: Convert PWM (Position Weight Matrix) to IUPAc sequences (java)
- ProcessedGOFormatter: Process SeqEnrich GO term output into a more usable format, filters p-values (old script, could do much better in R)

Methylation Project:
- 2017_08_01_bsmap_methylKit.R: Analyse data from BSMAP methylation call files with methylKit (Also create BED4 DMR files)
- bedtools_cmds.txt: Bash commands to create gene flanks and intersect gene/flanks with DMRs (Differentially Methylated Regions)
- extract_filter_intersectBed.R: Process the output of bedtools_cmds into a more usable format. Includes cmds for SeqEnrich prep
- methylation_heatmap_and_bargraph.R: Visualize methylation data by context, sample, and chromosome, visualize transposable element, gene, and CpG island density. (Mostly bar charts currently, heatmaps are not turning out great)

Pipeline: 
- SeqAnalysisPipeline: (perl) Guides user through the Belmonte lab's Sequencing Analysis Pipeline (Trimmomatic, TopHat2 with Bowtie2-build, Cuffquant, Cuffnorm, and CuffDiff) using cmd line prompts. Does preliminary error checking, saves a copy of the commands that were run, and organizes the output. Intructions/tutorial PDF file included in repository.
