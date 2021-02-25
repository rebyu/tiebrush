# Simulated Example Data
This directory contains a small simulated dataset for testing and getting familiar with TieBrush.

## Running TieBrush and TieCov

### Collapsing reads and building representative alignment for each tissue with TieBrush

```tiebrush -o t1/t1.bam t1/t1s0.bam t1/t1s1.bam t1/t1s2.bam t1/t1s3.bam t1/t1s4.bam t1/t1s5.bam t1/t1s6.bam t1/t1s7.bam t1/t1s8.bam t1/t1s9.bam```
```tiebrush -o t2/t2.bam t2/t2s0.bam t2/t2s1.bam t2/t2s2.bam t2/t2s3.bam t2/t2s4.bam t2/t2s5.bam t2/t2s6.bam t2/t2s7.bam t2/t2s8.bam t2/t2s9.bam```

### Extracting BED-formatted summaries for each tissue with TieCov

```tiecov -s t1/t1.sample.bed -c t1/tb.coverage.bed -j t1/tb.junctions.bed t1/tb.bam```
```tiecov -s t2/t2.sample.bed -c t2/tb.coverage.bed -j t2/tb.junctions.bed t2/tb.bam```

### Visualizing results in IGV
Upon completion, tiecov will have generated several summary files for each tissue.
Along with the collapsed alignment, these summaries can be viewed and scrutinized using
genome browsers such as [IGV](http://software.broadinstitute.org/software/igv/)
that support BAM (for collapsed alignment), BED (for the junctions track)
and BEDgraph (for the sample and coverage tracks) formats.

## Data
Example dataset contains mapped simulated reads (BAM) from two loci (NEFL and SLC25A3)
which are known to exhibit differential expression and splicing in several tissues.
Data provided is separated into two tissues (t1 and t2), each comprised of 10 samples.

NEFL (chr8:24,950,955-24,956,721 (GRCh38/hg38)) neurofilament-coding gene is known to be primarily expressed in certain cells of the 
Brain, however low expression are typically observed in several cell types in the heart.

![](https://github.com/alevar/tiebrush/blob/master/example/nefl.sim.png)
**Figure 1.** *Comparison of transcription of the NEFL gene in brain (top) and 
heart (bottom) tissues from the simulated dataset. Each tissue is represented by 
three tracks produced by TieCov: read coverage (top), percent of samples containing the 
reads (middle), and splice junctions (bottom). The plot illustrates the higher 
prevalence and expression of the gene in brain tissue*

SLC25A3 (chr12:98,593,591-98,606,379 (GRCh38/hg38)) can be expressed by several dominant
alternative isoforms, with splicing patterns being known to be tissue specific. In case of 
cells from heart and brain - alternative forms used differ in the third exon.

![](https://github.com/alevar/tiebrush/blob/master/example/slc25a3.sim.png)
**Figure 2.** *Comparison of transcription of the SLC25A3 gene in brain (top) 
and heart (bottom) tissues from the simulated dataset. While the gene is expressed 
in both tissues, coverage data clearly indicate an exon switch where the 3rd and 
4th exons are expressed at dramatically different levels in the two tissues*

## Simulation Protocol
1. 10 samples were randomly selected from brain tissue and another 10 samples were randomly selected 
from the heart tissue of the GTEx project (Table 1). 
2. Samples were aligned with HISAT2 against the GRCh.38(patch 12).
3. Transcriptomes for each sample were individually assembled and quantified using StringTie2.
4. Using GffCompare samples were individually compared against the RefSeq annotation
5. Transcript structures (GTF) and coverage data of isoforms overlapping the NEFL and SLC25A3 genes 
   were extracted for each sample.
6. Extracted GTF and coverage data were used to simulate reads using polyester.
7. Fasta files produced by polyester were directly converted into SAM format using sim2sam utility.

| Brain Samples | Heart Samples |
|---------------|:--------------|
|SRR1368222     | SRR1345221    |
|SRR1370173  	| SRR1441464    |
|SRR1374221 	| SRR1432366    |
|SRR1363810     | SRR661229     |
|SRR1467882     | SRR1337388    |
|SRR1471773  	| SRR1458484    |
|SRR1473590   	| SRR816358     |
|SRR814989   	| SRR1352213    |
|SRR817658   	| SRR1435584    |
|SRR818146   	| SRR1345736    |
**Table 1.** *GTEx samples used to model parameters for the simulation*