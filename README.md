# complex-organisms-WGS-pipeline
Reproducible GATK4 WGS variant-calling pipeline: SRA/FASTQ to filtered SNPs and indels. This pipeline employs fastp for QC and trimming, BWA-MEM2 for alignment, and GATK HaplotypeCaller with configurable ploidy for humans, animals and polyploid plants.
