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
* Variant filtering with `vcfFilter.pl` and `filterSomeMore.pl`

Used the following filters: 2X coverage (2496 reads), 10 alt. reads, not fixed, BQRS max abs. = 3, MQRS max abs. 2.5, RPRS max abs. 2, minimum ratio of varriant confidence to non-reference read depth (QD) 2, minimum mapping quality 30, missing data for fewer than 250 (80% with data), biallelic SNPs only.

Ended up with 163,850 SNPs in /uufs/chpc.utah.edu/common/home/gompert-group2/data/dimension_med_gbs/Variants/filtered2x_msativa_cg.vcf. There are also files for individual chromosomes (lg*.vcf).

Next, I dropped SNPs with > mean + 3 SD coverage (possible repeats). This left 161,008 SNPs in morefilter_filtered2x_msativa_cg.vcf.

2. *L. melissa*

* Filtered for PhiX, Split files, then demultiplexed as per usual. PhiX reference is /uufs/chpc.utah.edu/common/home/u6000989/data/phixSeqIndex/NC_001422.fasta.

* Indexed new version of *L. melissa* genome, i.e., dovetail (Hi-C plus Chicago) with PacBio gap filling: /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta

````bash
bwa index -a bwtsw Lmel_dovetailPacBio_genome.fasta
````
 * DNA sequences aligned to the *L. melissa* reference genome (see above) with `bwa` (0.7.17-r1188) with `bwa mem`

```bash
bwa mem -t 1 -k 15 -r 1.3 -T 30 -R '@RG\tID:Lmelissa-ID\tPL:ILLUMINA\tLB:Lmelissa-ID\tSM:L.melissa-ID' /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta > alnID.sam
```
 * Compress, sort and index with `samtools` (version 1.10)

```bash
samtools view -b -O BAM -o ID.bam ID.sam
samtools sort -O BAM -o ID.sorted.bam ID.bam
samtools index -b ID.sorted.bam
```

* Variant calling with `bcftools` version 1.9

```bash
bcftools mpileup -C 50 -d 250 -f /uufs/chpc.utah.edu/common/home/gompert-group3/data/LmelGenome/Lmel_dovetailPacBio_genome.fasta -q 30 -Q 20 -I -b lmel_bams.txt -o lmel_variants.bcf -O u -a FORMAT/AD,FORMAT/DP
bcftools call -v -c -p 0.01 -P 0.001 -O v -o lmel_variants.vcf lmel_variants.bcf 
```

* Variant filtering with `vcfFilter.pl` and `filterSomeMore.pl`

Used the following filters: 2X coverage (2302 reads), 10 alt. reads, not fixed, Man-Whitney P for BQB = 0.01, Man-Whitney P for RPB = 0.01, minimum mapping quality 30, missing data for fewer than 230 (80% with data), biallelic SNPs only.

Ended up with 64,061 SNPs in /uufs/chpc.utah.edu/common/home/gompert-group2/data/dimension_lyc_gbs/Variants/filtered2x_lmel_variants.vcf. 

Next, I dropped SNPs with > mean + 3 SD coverage (possible repeats). This left 63,194 SNPs in morefilter_filtered2x_lmel_variants.vcf.

# Estimating genotypes

1. *M. sativa*

* Use LDA of PCA to generate initial values for `entropy` (version 2.0)

````R
library(data.table)
g<-as.matrix(fread("pntest_filtered2x_msativa.txt",header=FALSE))

## pca on the genotype covariance matrix
pcgcov<-prcomp(x=t(g),center=TRUE,scale=FALSE)

## kmeans and lda
library(MASS)
k2<-kmeans(pcgcov$x[,1:5],2,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k3<-kmeans(pcgcov$x[,1:5],3,iter.max=10,nstart=10,algorithm="Hartigan-Wong")

ldak2<-lda(x=pcgcov$x[,1:5],grouping=k2$cluster,CV=TRUE)
ldak3<-lda(x=pcgcov$x[,1:5],grouping=k3$cluster,CV=TRUE)

write.table(round(ldak2$posterior,5),file="ldak2.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak3$posterior,5),file="ldak3.txt",quote=F,row.names=F,col.names=F)
````
* Estimate genotypes with `entropy`, with k = 2,3 and 12 chains each. Here is one example:

````bash
entropy_mp ../Variants/filtered2x_msativa.gl -n 4 -m 1 -l 1000 -b 500 -t 5 -k 2 -out_msat_K2_ch0.hdf5  -q ldak2.txt -s 20
````
* Then summarize the posterior with `estpost.entropy` version 2.0.

```bash
estpost.entropy -p gprob -s 0 -w 0 -o msat_geno.txt out_msat_k*hdf5
# parameter dimensions for gprob: loci = 161008, ind = 1248, genotypes = 5, chains = 24
```

2. *L. melissa*

* Use LDA of PCA to generate initial values for `entropy` (version 2.0)

````R
## see initq_lmel.R
library(data.table)
g<-as.matrix(fread("pntest_melissa_filtered2x_lmelissa.txt",header=FALSE))

## pca on the genotype covariance matrix
pcgcov<-prcomp(x=t(g),center=TRUE,scale=FALSE)

## kmeans and lda
library(MASS)
k2<-kmeans(pcgcov$x[,1:5],2,iter.max=10,nstart=10,algorithm="Hartigan-Wong")
k3<-kmeans(pcgcov$x[,1:5],3,iter.max=10,nstart=10,algorithm="Hartigan-Wong")

ldak2<-lda(x=pcgcov$x[,1:5],grouping=k2$cluster,CV=TRUE)
ldak3<-lda(x=pcgcov$x[,1:5],grouping=k3$cluster,CV=TRUE)

write.table(round(ldak2$posterior,5),file="ldak2_lmel.txt",quote=F,row.names=F,col.names=F)
write.table(round(ldak3$posterior,5),file="ldak3_lmel.txt",quote=F,row.names=F,col.names=F)

save(list=ls(),file="initq_lmel.rdat")
````
* Estimate genotypes with `entropy`, with k = 2,3 and 12 chains. Here is one example:

````bash
entropy_mp ../Variants/filtered2x_msativa.gl -n 2 -m 1 -l 2000 -b 1000 -t 5 -k 2 -out_K2_ch0.hdf5  -q ldak2_lmel.txt -s 20
````
* Then summarize the posterior with `estpost.entropy` version 2.0.

```bash
estpost.entropy -p gprob -s 0 -w 0 out_k*hdf5
```

# Desribing genetic variation

I summarized patterns of genetic variation for the *M. sativa* and *L. melissa* individuals for the experiment based on the genotype data. This was done in `R`, see [plotPCAwMap.R](plotPCAwMap.R). PCA was perofrmed on the centered (but not standaridzed) gentoypes. Point estimates of allele frequencies for each source population were calcualted as the mean genotype for each locus and population divided by two. These were then used to estimate Nei's Gst, i.e.,:

```{R}
## fst (Nei Gst)
fst<-function(P=NA){
        ## P = allele freq. matrix, cols=pops, rows=SNPs
        H<-P * (1-P) * 2
        Hw<-apply(H,1,mean)## mean across pops for each locus
        pbar<-apply(P,1,mean)
        Ht<-pbar * (1-pbar) * 2
        Fst<-1 - mean(Hw)/mean(Ht)
        return(Fst)
}
```
I also computed LD between pairs of linked SNPs and summarized this at different distances. LD drops off quickly in *L. melissa* and is low in general; LD is a bit higher in *M. sativa* but still decays quickly to background levels [caldLD.R](caldLD.R).

# Preparing the phenotypic data

* Phenotypic data are in /uufs/chpc.utah.edu/common/home/gompert-group2/projects/dimensions_cg_experiment/Pheno/. These include chemistry data (log2, 1750 features plus 6 PCs), ...

* The main analyses use the data from the field plot. Data from the science garden and BST, which were used to validate/test our models, will be dealt with below. I used Moran's eigenvector maps to remove possible effects of space (location) within the common gardent. This essentially involves creating spatial variables based on a PCoA of a truncated (at nearest neighbors) Euclidean distane matrix, see for example [Dray et al., 2006](https://www.sciencedirect.com/science/article/pii/S0304380006000925?casa_token=DIpYsGyIWGAAAAAA:vJCYcARQofbW-Q_GZ3koZpW47d3kj88TlmVKg8nB5p-eysSmKrXbG7aw_MDc-qY5fkkkHsfNvg). I then used forward selection of variables following [Blanchet et al., 2008](https://esajournals.onlinelibrary.wiley.com/doi/pdf/10.1890/07-0986.1?casa_token=EC3ZTIn92BoAAAAA:B3sPgrMDh5WOBS8s6mecbJF34wl1jZseYoRa87orTzmzFceYwdmeNnX573UjI1RnlQ9ixLwztF1eWQ) to select spatial variables (eigenvectors) that explained each trait. Scaled residulas from these analyses (along with scaled versions of the original variables) were then used for analysis. See [summarizeFormatChem.R](summarizeFormatChem.R).

* Complementary analyses/code for cleaning and formatting caterpillar perfomance data, including removing the effect of hatch date, are [cleanCats.R](cleanCats.R) and  [summarizeFormatCats.R](summarizeFormatCats.R).

# Estimating PVE and polygenic scores for caterpillar performance with using BSLMMs with gemma

* I used `gemma` (version 0.95a) to estimate PVE and polygenic scores for the caterpillar performance traits based on *M. sativa* genotypes, *L. melissa* genotypes and genetic data from both species combined. Here is an example perl fork script for this that was used to fit models for the actual data set and a randomized version of the data set. The basic parameters, burnin = , length of chain = and running 10 chains total were used in all cases.

```{perl}
#!/usr/bin/perl
#
# fit gemma BSLMM for M. sativa, non-chemistry, real and randomized 
#
use Parallel::ForkManager;
my $max = 80;
my $pm = Parallel::ForkManager->new($max);

$g = "msat_geno";

foreach $p (@ARGV){
        open(WC, "head -n 1 $p | wc|");
        $wc = <WC>;
        @wc = split(/\s+/,$wc);
        $Nph = $wc[2];
        $p =~ m/^([a-zA-Z_]+)/;
        $base = $1;

        foreach $ph (1..$Nph){ 
                foreach $ch (0..9){
                        sleep 2;
                        $pm->start and next;
                        $o = "o_msat_fit_$base"."_ph$ph"."_ch$ch";
                        system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -maf 0 -w 200000 -s 1000000\n";
                        $pm->finish;
                }
        }
}
$pm->wait_all_children;
```

* In each case, I then summarized hyperparameter posteriors and obtained estimates of model-averaged effects. See the example below, and [calpost.pl](calpost.pl) and [grabMavEffects.pl](grabMavEffects.pl)

```bash
perl calpost.pl o_lmel_fit_gemmalmel_pheno_residTraits_ph*ch0.hyp.txt
perl calpost.pl o_lmel_fit_gemmalmel_pheno_rawTraits_ph*ch0.hyp.txt
perl calpost.pl o_lmel_fit_randomlmel_pheno_rawTraits_ph*ch0.hyp.txt
perl calpost.pl o_lmel_fit_randomlmel_pheno_residTraits_ph*ch0.hyp.txt

perl grabMavEffects.pl o_lmel_fit_randomlmel_pheno_residTraits_ph*ch0.param.txt
```
* Next, I estimated polygenic scores for the performance traits based on plant, caterpillar or combined genetics. I did this using model-averaged effect estimates (i.e., all genetic loci contributed to some extent). For details of this calculation, along with some summary plots, see [mkGenArchPlots.R](mkGenArchPlots.R). 

# Within garden genomic prediction and cross-validation

I used 10-fold cross-validation to predict caterpillar performance from *M. sativa* genetics, *L. melissa* genetics, and *M. sativa* + *L. melissa* genetics. This first involved creating 10 data sets for each variable each with a random ~10% of individuals with data set to NA (the test set; the other 90% served as the training set). See [mkCvPheno.R](mkCvPheno.R).

I then ran `gemma` 10 times per original performance trait, that is, once per cross-validation data set with only a single chain in each case. This included the model fit and genomic prediction for the missing (test) individuals. See, e.g., the code blow for the case of *M. sativa* (scripts for other genetic data were essentially identical). 

```{perl}
use Parallel::ForkManager;
my $max = 80;
my $pm = Parallel::ForkManager->new($max);

$g = "msat_geno";

foreach $p (@ARGV){
	open(WC, "head -n 1 $p | wc|");
	$wc = <WC>;
	@wc = split(/\s+/,$wc);
	$Nph = $wc[2];
	$p =~ m/^([a-zA-Z_]+)/;
	$base = $1;

	foreach $ph (1..$Nph){ 
		$ch = 0;
		sleep 2;
		$pm->start and next;
		$o = "o_msat_fit_$base"."_ph$ph"."_ch$ch";
    		system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -outdir output_cv -maf 0 -w 200000 -s 1000000\n";
    		system "gemma -g $g -p $p -n $ph -o $o -maf 0 -epm output_cv/$o.param.txt -emu output_cv/$o.log.txt -predict 1 -o pred_$o -outdir output_cv \n";
		$pm->finish;
	}
}
$pm->wait_all_children;
```

# Predictive power between gardens

I next tested our ability to predict caterpillar performance trait values from genotype by generating genomic predictions of performance for caterpillars reared on *M. sativa* from a second, smaller common garden comprising 180 *M. sativa*. This second garden, planted on the USU campus (i.e., the [Science Garden](https://biology.usu.edu/sciencegarden/)), included plants from six of the 11 *M. sativa* source sites (30 plants per site from six maternal families) and caterpillars from each of the sites used in the main experiment (ALP, APLL, AWFS, BST, VUH, and VIC).

I calculated genomic predictions of performance trais based on from model-avereage effects for plant, caterpillar or plant and caterpillar genetics inferred from the main garden and genotype estimates for the 172 plants and 156 caterpillars successfully sequenced from the common gardent. I then computed the Pearson correlation between observed (after removing effects of hatch date and block) and predicted performance for each trait/genetic data set as a metric of predictive perofrmance. See [GenPredsSG.R](GenPredsSG.R).

# USU greenhouse experiment

We conducted another complimentary experiment using greenhouse-grown *M. sativa* in 2018. This experiment again used plants from ALP, APLL, AWFS, BST, VUH and VIC (1001 total), but included caterpillars from four *L. melissa* populations (6-13 females per site, 672 caterpillars total), *Colias eurytheme* from the Greenville Farm garden (20 females caught, 133 caterpillars total) and *Vanessa cardui* from Carolina Biological (five units of 30-35 eggs, 196 caterpillars total). Lauren's initial notes on this experiment are on the [old site](https://sites.google.com/site/gompertlabnotes/home/researcher-pages/lauren-lucas/alfalfa-2018). Tara completed many initial analyses of this experiment as well. Here, I am focused on two questions: (i) does among-population plant (and caterpillar for *L. melissa*) variation matter for 8 and 14 day weight (performance), and (ii) does plant genotype (population) have consistent effects of caterpillar performance across different butterfly populations and species. I addressed (i) by estimating variance components for plant (and caterpillar) population and (ii) by computing correlations in plant effects across butterfly populations and species.  See [VarComps.R](VarComps.R) for details, but the short answer is that plant genetics matters (more for 14 day weight), caterpillar genetics matters for 8 day weight in *L. melissa*, and plant population effects on performance are remarkably consitent across butterfly populations and species.

# Plant trait mapping

I used the same BSLMM approach to map 1760 plant traits from the Greenville Experimental Farm garden; *N* = 1080 individuals. These included 9 morphological traits: leaf length, leaf width, leaf area, leaf shape, leaf dry weight, SLA, trichome density, leaf toughness and plant height ; field herbivory; and 1750 chemical features (led by Casey at UNR). Mearements data are all from summer 2019. See [notes and data from google drive](https://drive.google.com/drive/folders/1gEP_t7SxYBe579Sj8JPUGOea6yxtyHQS?usp=sharing). Data are in `/uufs/chpc.utah.edu/common/home/gompert-group2/projects/dimensions_cg_experiment/Pheno`.

First, I removed spatial effects on plant traits using stepwise regression with Moran's eigenvector spatial coordinates as was done for caterpillar performance. See [summarizeFormatChem.R](summarizeFormatChem.R) and [summarizeFormatPlantTr.R](summarizeFormatPlantTr.R). Then mapping was done with `gemma` version (). This always invovled 10 chains each of 1 million steps and a 200,000 step burnin. Here is an example of the `perl` conrol script (done in many batches as this was huge):

```{perl}
#!/usr/bin/perl

use Parallel::ForkManager;
#my $max = 50;
my $max = 60;
my $pm = Parallel::ForkManager->new($max);

$g = "msat_geno";
$ph_lb = 1701;
$ph_ub = 1750;
#$dir = "output_rawchem_$ph_ub";
$dir = "output_ranchem_$ph_ub";
#$dir = "output_chem_$ph_ub";

foreach $p (@ARGV){
        $p =~ m/^([a-zA-Z_]+)/;
        $base = $1;

        foreach $ph ($ph_lb .. $ph_ub){ 
                foreach $ch (0..9){
                        sleep 2;
                        $pm->start and next;
                        $o = "o_msat_fit_$base"."_ph$ph"."_ch$ch";
                        system "gemma -g $g -p $p -bslmm 1 -n $ph -o $o -maf 0 -w 200000 -s 1000000 -outdir $dir\n";
                        $pm->finish;
                }
        }
}
```
This entire procedure was repeated for 1760 randomized data sets, one randomized data set for each of the original 1760 traits. Bulk MCMC results are in teh `output_chem*` and `output_ranchem` subdirectories. Posteriores were summarized exactly as was done for the performance traits. Model-averaged estimates of of SNP effects for all traits are in `/uufs/chpc.utah.edu/common/home/gompert-group2/projects/dimensions_cg_experiment/Gemma/files_mav`. These were used to compute polygenic scores, see [computePolyScores.R](computePolyScores.R).

# LASSO regression models for caterpillar performance as a function of plant-trait polygenic scores

We next used LASSO regression to fit models for the polygenic scores for each caterpillar performance trait (based on *M. sativa*) genetics as a function of the polygenic scores from all 1760 plant traits. This was done with `glmnet` (version  4.0-2) in `R` and used 10-fold cross-validation to selection lambda. I also measured both variance explained by the model and the 10-fold cross-validation predictive power. See [glmnetPerformPS.R](glmnetPerformPS.R). Variance explained was compared for the observed versus randomized response variables.

As a test of how many *independent* genetic factors (from *M. sativa*) contributed to the caterpillar performance polygenic scores, I also fit LASSO models on PCs of the plant trait and chemistry polygenic scores. These gave generally similar results to the analyses above providing (I think) pretty clear evidence of many genetic variants involved in different traits contributing [glmnetPerformPS_PCA.R](glmnetPerformPS_PCA.R).

I repeated this analysis with observed performance values (residuals after hatch data and space, but not polygenic scores). Here too randomized data sets were also analyzed. See [glmnetPSObs.R](glmnetPSObs.R).

Next, I fit LASSO regression models that allowed for interactions between plant-trait polygenic scores and caterpillar genotype. For *L. melissa* caterpillar genotype I used the first four PCs from a PCA of the caterpillar genotype matrix. These PCs were also included as main effects. See [glmnetPerformPS_epistasis.R](glmnetPerformPS_epistasis.R). Finally, I repeated these analyses with the observed performance data rather than polygenic scores (but still using residuals) to make sure the poor fit of the interaction models was not simply due to the additive models used to estimate performance polygenic scores. This gave generally similar results, with not overtall increase in explanatory power for models with epsistasis. See [glmnetPerformObs_epistasis.R](glmnetPerformObs_epistasis.R).
