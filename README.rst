TieBrush and TieCov - efficient methods for aggregating and summarizing aligned sequences across large datasets
===============================================================================================================

.. image:: https://img.shields.io/badge/License-MIT-blue.svg
    :target: https://opensource.org/licenses/MIT
    :alt: MIT License

.. contents::
   :local:
   :depth: 2

Introduction
^^^^^^^^^^^^

TieBrush is a simple yet efficient method for merging redundant information from multiple alignment files, 
designed to enable rapid manipulation of extremely large sequencing datasets. The method is specifically 
designed to optimize investigations of RNA, whole-genome, exome and other types of sequencing experiments. 
TieBrush preserves much of the original information in a greatly condensed representation as a BAM file, 
which allows manipulation and extraction of dataset and subset-specific statistics using tools within 
the package, and which is also compatible with other common utilities.

This utility aims to merge/collapse "duplicate" read alignments (same location with the same CIGAR string),
across multiple sequencing samples (multiple input BAM files), adding custom SAM tags in order to keep
track of the "alignment multiplicity" count (how many times the same alignment is seen across all
input data) and "sample count" (how many samples show that same alignment).
The initial goal is to generate this composite BAM file which multiplexes read alignments
from many sequencing samples, painting a comprehensive "background" picture of read alignments
with their counts across many samples.

Installation
^^^^^^^^^^^^

Building from source
""""""""""""""""""""

If you want to build it from source, we recommend cloning the git repository as shown below to ensure
fixed releases of any dependencies are fetched and compiled with the software.

::

    $ git clone https://github.com/alevar/tiebrush.git --recursive
    $ cmake -DCMAKE_BUILD_TYPE=Release .
    $ make -j4
    $ sudo make install

If you are using a very old version of Git (< 1.6.5) the flag ``--recursive`` does not exist.
In this case you need to clone the submodule separately  (``git submodule update --init --recursive``).

**Requirements**

Operating System
  GNU/Linux, Mac

Compiler
  GCC ≥ 4.9, LLVM/Clang ≥ 3.8

Build system
  CMake ≥ 2.8

Language support
  C++11

Methods
^^^^^^^

TieBrush
""""""""

Summarize and filter read alignments from multiple sequencing samples (taken as sorted BAM files).
This utility aims to merge/collapse "duplicate" read alignments (same location with the same
CIGAR string), across multiple sequencing samples (multiple input BAM files), adding custom SAM tags
in order to keep track of the "alignment multiplicity" count (how many times the same alignment is
seen across all input data) and "sample count" (how many samples show that same alignment).

The goal is to generate this composite BAM file which multiplexes read alignments from many sequencing
samples, painting a comprehensive "background" picture of read alignments with their counts across
many samples.

  tiebrush [-o <outfile>.bam] list.txt | in1.bam in2.bam ... inN.bam

  -C, --cigar        merge if only CIGAR string is the same.
  -P, --clip         merge if clipped CIGAR string is the same.
  -E, --exon         merge if exon boundaries are the same.
  -S, --keep_supp    keep supplementary alignments.
  -M, --keep_unmap   keep unmapped reads as well.
  -K, --keep_names   keep names in the original format. if not set - all values will be set incrementally.
  -U, --keep_quals   keep quality strings for the collapsed records. Note that quality strings are randomly chosen if 2 or more records are collapsed.
  -N                maximum NH score (if available) to include.
  -Q                minimum mapping quality to include.
  -F                bits in SAM flag to use in read comparison.

SAM tags implemented
--------------------
 1. __YC__:i:N stores the number of alignments that were merged into this alignment record (multiplicity count)
 2. __YX__:i:N stores the number of samples that have this alignment (sample count)
 3. __YD__:i:N keeps track of the maximum number of contiguous bases preceding the start of the read alignment in the samples(s) that it belongs to. In other words, if the current alignment is part of an exon-overlapping bundle (strand specific!), this value holds the maximum distance from the beginning of the bundle to the start of this alignment, across all samples having this alignment. If the alignment is not in a bundle (i.e. it is preceded by a uncovered region as it is not overlapped by any another alignment with a lower start position), in all the individual samples where that alignment is present, then the YD value is 0 and the tag is omitted from the output file produced by TieBrush. That means that all the alignments lacking a YD tag in the TieBrush output start at the very beginning of an exon-overlapping bundle (i.e. are not overlapped by a preceding alignment with a lower start coordinate).

If either YC or YX tags are missing (i.e. GBamRecord::__tag_int__() call returns 0) then the alignment is unique (when YC is 0) or only one sample has it (if YX is 0). The actual count in these cases is obviously 1.

TieCov
""""""

The tiecov utility can take the output file produced by TieBrush and can generate the following auxiliary base/junction coverage files:
 1. a BedGraph file with the coverage data (see http://genome.ucsc.edu/goldenPath/help/bedgraph.html); this file can be converted to BigWig (using bedGraphToBigWig) or to TDF format (using igvtools) in order to be loaded in IGV as an additional coverage track
 2. a junction BED file which can be loaded directly in IGV as an additional junction track (http://software.broadinstitute.org/software/igv/splice_junctions)
 3. a heatmap BED that uses color intensity to represent the number of samples that contain each position.

  tiecov [-b out.flt.bam] [-s out.sample.bed] [-c out.coverage.bedgraph] [-j out.junctions.bed] in.bam

  -b    bam file after applying filters (-N/-Q)
  -s    BED file with number of samples which contain alignments for each interval.
  -c    BedGraph file with coverage for all mapped bases.
  -j    BED file with coverage of all splice-junctions in the input file.
  -N    maximum NH score (if available) to include when reporting coverage
  -Q    minimum mapping quality to include when reporting coverage

TieWrap
"""""""

TieWrap is a small utility script provided to make running tiebrush on large datasets a bit easier.
Unlike TieBrush, TieWrap can be launched with as many input files as needed and will automatically
divide them into batches processing and combining batches to produce a single representation at the end.
All standard TieBrush arguments can be passed over to TieWrap. Additionally size of individual batches
as well as the concurrency parameters can be set explicitely.

  tiewrap.py [-h] -o OUTPUT [-C] [-P] [-E] [-S] [-M] [-N MAX_NH] [-Q MAX_MAP_QUAL] [-F FLAGS] [-t THREADS] [-b BATCH_SIZE] list.txt | in1.bam in2.bam ... inN.bam

  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        output file
  -C, --cigar           merge if only CIGAR string is the same
  -P, --clip            merge if clipped CIGAR string is the same
  -E, --exon            merge if exon boundaries are the same
  -S, --keep-supp       keep supplementary alignments
  -M, --keep-unmap      keep unmapped reads
  -N MAX_NH, --max-nh MAX_NH
                        maximum NH score of the reads to retain
  -Q MAX_MAP_QUAL, --max-map-qual MAX_MAP_QUAL
                        maximum NH score of the reads to retain
  -F FLAGS, --flags FLAGS
                        bits in SAM flag to use in read comparison
  -t THREADS, --threads THREADS
                        number of threads to use
  -b BATCH_SIZE, --batch-size BATCH_SIZE
                        Number of input files to process in a batch on each
                        thread
