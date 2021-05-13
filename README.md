# DimensionsExperiment
Code and overview for DoB common garden experiment

# Raw data

DNA sequences for *L. melissa*: /uufs/chpc.utah.edu/common/home/gompert-group2/data/dimension_lyc_gbs

This includes plates Lyc-048 to Lyc-059 (libraries TS-1 to TS-3).

DNA sequences for *M. sativa*: /uufs/chpc.utah.edu/common/home/gompert-group2/data/dimension_med_gbs

This includes plates Alf-001 to Alf-014 (libraries TS-3 to TS-7).

Plant chemistry data: /uufs/chpc.utah.edu/common/home/gompert-group2/data/dimensions_chem

# DNA sequence filtering, alignment, variant calling and filtering

1. *M. sativa*

  * DNA sequences aligned to the *M. sativa* reference genome with *bwa* (0.7.17-r1188) with *bwa mem*

```bash
bwa mem -t 1 -k 15 -r 1.3 -T 30 -R '@RG\tID:Msativa-ID\tPL:ILLUMINA\tLB:Msativa-ID\tSM:Msativa-ID' /uufs/chpc.utah.edu/common/home/gompert-group2/data/alfalfa_genome/sc_final_genome.fasta ID.fastq > alnID.sam
```
