# DimensionsExperiment
Code and overview for DoB common garden experiment

# Raw data

DNA sequences for *L. melissa*: /uufs/chpc.utah.edu/common/home/gompert-group2/data/dimension_lyc_gbs

This includes plates Lyc-048 to Lyc-059 (libraries TS-1 to TS-3).

DNA sequences for *M. sativa*: /uufs/chpc.utah.edu/common/home/gompert-group2/data/dimension_med_gbs

This includes plates Alf-001 to Alf-014 (libraries TS-4 to TS-7).

Plant chemistry data: /uufs/chpc.utah.edu/common/home/gompert-group2/data/dimensions_chem

# DNA sequence filtering, alignment, variant calling and filtering

1. *M. sativa*

 * Filtered for PhiX, Split files, then demultiplexed as per usual. PhiX reference is /uufs/chpc.utah.edu/common/home/u6000989/data/phixSeqIndex/NC_001422.fasta.

 * DNA sequences aligned to the *M. sativa* reference genome with `bwa` (0.7.17-r1188) with `bwa mem`

```bash
bwa mem -t 1 -k 15 -r 1.3 -T 30 -R '@RG\tID:Msativa-ID\tPL:ILLUMINA\tLB:Msativa-ID\tSM:Msativa-ID' /uufs/chpc.utah.edu/common/home/gompert-group2/data/alfalfa_genome/sc_final_genome.fasta ID.fastq > alnID.sam
```

 * Compress, sort and index with `samtools` (version 1.10)

```bash
samtools view -b -O BAM -o ID.bam ID.sam
samtools sort -O BAM -o ID.sorted.bam ID.bam
samtools index -b ID.sorted.bam
```
* Variant calling with `GATK` (version 4.1) (necessary because *M. sativa* is tetraploid)

```bash
## generate g.vcf files (rep. across individuals)
gatk --java-options "-Xmx20g" HaplotypeCaller -R uufs/chpc.utah.edu/common/home/gompert-group2/data/alfalfa_genome/sc_final_genome.fasta  -I PATH/ID.sorted.bam -O PATH/ID.sorted.g.vcf -heterozygosity 0.001 -mbq 20 -ERC GVCF -ploidy 4
## combine g.vcf files (two steps, one example of first shown)
gatk --java-options "-Xmx40g" CombineGVCFs -R /uufs/chpc.utah.edu/common/home/gompert-group2/data/alfalfa_genome/sc_final_genome.fasta --variant aln_FP-10-AFAL-10-1.sorted.g.vcf --variant aln_FP-10-AFAL-12-1.sorted.g.vcf --variant aln_FP-10-AFAL-15-2.sorted.g.vcf --variant aln_FP-10-AFAL-17-1.sorted.g.vcf --variant aln_FP-10-AFAL-18-4.sorted.g.vcf --variant aln_FP-10-AFAL-6-3.sorted.g.vcf --variant aln_FP-10-ALP-10-1.sorted.g.vcf --variant aln_FP-10-ALP-14-1.sorted.g.vcf --variant aln_FP-10-ALP-15-5.sorted.g.vcf --variant aln_FP-10-ALP-18-5.sorted.g.vcf --variant aln_FP-10-ALP-19-2.sorted.g.vcf --variant aln_FP-10-ALP-19-3.sorted.g.vcf --variant aln_FP-10-ALP-20-4.sorted.g.vcf --variant aln_FP-10-ALP-21-4.sorted.g.vcf --variant aln_FP-10-ALP-6-2.sorted.g.vcf --variant aln_FP-10-ALP-6-3.sorted.g.vcf --variant aln_FP-10-ALP-6-5.sorted.g.vcf --variant aln_FP-10-APLL-10-1.sorted.g.vcf --variant aln_FP-10-APLL-10-4.sorted.g.vcf --variant aln_FP-10-APLL-17-3.sorted.g.vcf --variant aln_FP-10-APLL-2-1.sorted.g.vcf --variant aln_FP-10-APLL-3-3.sorted.g.vcf --variant aln_FP-10-APLL-5-2.sorted.g.vcf --variant aln_FP-10-APLL-9-3.sorted.g.vcf --variant aln_FP-10-AWFS-15-1.sorted.g.vcf --variant aln_FP-10-AWFS-17-2.sorted.g.vcf --variant aln_FP-10-AWFS-3-1.sorted.g.vcf --variant aln_FP-10-AWFS-9-2.sorted.g.vcf --variant aln_FP-10-BST-17-2.sorted.g.vcf --variant aln_FP-10-BST-18-4.sorted.g.vcf --variant aln_FP-10-BST-28-3.sorted.g.vcf --variant aln_FP-10-BST-28-5.sorted.g.vcf --variant aln_FP-10-CKV-01-3.sorted.g.vcf --variant aln_FP-10-CKV-02-5.sorted.g.vcf --variant aln_FP-10-CKV-06-1.sorted.g.vcf --variant aln_FP-10-CKV-07-5.sorted.g.vcf --variant aln_FP-10-CKV-09-1.sorted.g.vcf --variant aln_FP-10-CKV-12-5.sorted.g.vcf --variant aln_FP-10-CKV-13-1.sorted.g.vcf --variant aln_FP-10-CKV-26-1.sorted.g.vcf --variant aln_FP-10-CKV-5-3.sorted.g.vcf --variant aln_FP-10-FRM-10-5.sorted.g.vcf --variant aln_FP-10-FRM-13-5.sorted.g.vcf --variant aln_FP-10-FRM-1-3.sorted.g.vcf --variant aln_FP-10-FRM-1-5.sorted.g.vcf --variant aln_FP-10-FRM-18-4.sorted.g.vcf --variant aln_FP-10-FRM-25-3.sorted.g.vcf --variant aln_FP-10-FRM-27-3.sorted.g.vcf --variant aln_FP-10-HBR-12-5.sorted.g.vcf --variant aln_FP-10-HBR-3-5.sorted.g.vcf --variant aln_FP-10-HBR-8-4.sorted.g.vcf --variant aln_FP-10-LIK-10-5.sorted.g.vcf --variant aln_FP-10-LIK-3-1.sorted.g.vcf --variant aln_FP-10-LIK-3-2.sorted.g.vcf --variant aln_FP-10-LIK-7-4.sorted.g.vcf --variant aln_FP-10-VIC-13-1.sorted.g.vcf --variant aln_FP-10-VIC-14-4.sorted.g.vcf --variant aln_FP-10-VIC-18-4.sorted.g.vcf --variant aln_FP-10-VIC-19-2.sorted.g.vcf --variant aln_FP-10-VIC-2-1.sorted.g.vcf --variant aln_FP-10-VIC-22-5.sorted.g.vcf --variant aln_FP-10-VIC-3-4.sorted.g.vcf --variant aln_FP-10-VIC-5-2.sorted.g.vcf --variant aln_FP-10-VIC-8-2.sorted.g.vcf --variant aln_FP-10-VUH-13-5.sorted.g.vcf --variant aln_FP-10-VUH-14-4.sorted.g.vcf --variant aln_FP-10-VUH-14-5.sorted.g.vcf --variant aln_FP-10-VUH-17-3.sorted.g.vcf --variant aln_FP-10-VUH-25-5.sorted.g.vcf --variant aln_FP-10-VUH2-7-4.sorted.g.vcf --variant aln_FP-10-VUH-5-2.sorted.g.vcf --variant aln_FP-11-AFAL-10-4.sorted.g.vcf --variant aln_FP-11-AFAL-14-5.sorted.g.vcf --variant aln_FP-11-AFAL-2-2.sorted.g.vcf --variant aln_FP-11-AFAL-5-4.sorted.g.vcf --variant aln_FP-11-AFAL-6-1.sorted.g.vcf --variant aln_FP-11-ALP-13-4.sorted.g.vcf --variant aln_FP-11-ALP-19-4.sorted.g.vcf --variant aln_FP-11-ALP-22-4.sorted.g.vcf --variant aln_FP-11-ALP-23-2.sorted.g.vcf --variant aln_FP-11-ALP-26-3.sorted.g.vcf --variant aln_FP-11-ALP-26-5.sorted.g.vcf --variant aln_FP-11-ALP-27-3.sorted.g.vcf --variant aln_FP-11-ALP-6-4.sorted.g.vcf --variant aln_FP-11-APLL-13-2.sorted.g.vcf --variant aln_FP-11-APLL-14-1.sorted.g.vcf --variant aln_FP-11-APLL-14-2.sorted.g.vcf --variant aln_FP-11-APLL-15-3.sorted.g.vcf --variant aln_FP-11-APLL-16-3.sorted.g.vcf --variant aln_FP-11-APLL-8-2.sorted.g.vcf --variant aln_FP-11-APLL-8-4.sorted.g.vcf --variant aln_FP-11-AWFS-10-5.sorted.g.vcf --variant aln_FP-11-AWFS-13-4.sorted.g.vcf --variant aln_FP-11-AWFS-14-4.sorted.g.vcf --variant aln_FP-11-AWFS-19-5.sorted.g.vcf --variant aln_FP-11-AWFS-22-1.sorted.g.vcf -O combinded_1.g.vcf
gatk --java-options "-Xmx540g" CombineGVCFs -R /uufs/chpc.utah.edu/common/home/gompert-group2/data/alfalfa_genome/sc_final_genome.fasta --variant combinded_1.g.vcf --variant combinded_2.g.vcf --variant combinded_3.g.vcf --variant combinded_4.g.vcf --variant combinded_5.g.vcf --variant combinded_6.g.vcf --variant combinded_7.g.vcf --variant combinded_8.g.vcf --variant combinded_9.g.vcf --variant combinded_10.g.vcf --variant combinded_11.g.vcf --variant combinded_12.g.vcf --variant combinded_13.g.vcf -O combinded_med_sativa.g.vcf
## joint variant calling (seperately for each chromosome, here is chrom. 1)
gatk --java-options "-Xmx48g" GenotypeGVCFs -R /uufs/chpc.utah.edu/common/home/gompert-group2/data/alfalfa_genome/sc_final_genome.fasta  --heterozygosity 0.001 --intervals 1 --V combinded_med_sativa.g.vcf -O combinded_med_sativa_lg1.vcf
## combined  to create combined_SG.g.vcf
```
* Variant filtering with `vcfFilter.pl`

Used the following filters: 2X coverage (2496 reads), 10 alt. reads, not fixed, BQRS max abs. = 3, MQRS max abs. 2.5, RPRS max abs. 2, minimum ratio of varriant confidence to non-reference read depth (QD) 2, minimum mapping quality 30, missing data for fewer than 250 (80% with data), biallelic SNPs only.

Ended up with 163,850 SNPs in /uufs/chpc.utah.edu/common/home/gompert-group2/data/dimension_med_gbs/Variants/filtered2x_msativa_cg.vcf. There are also files for individual chromosomes (lg*.vcf).

