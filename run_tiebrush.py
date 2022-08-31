#!/usr/bin/env python

#===================================================================
# run tiebrush given input file list
#    python run_tiebrush.py --input /path/to/input_csv
#                           --outdir /path/to/output/directory
#                           --threads 10 
#                           --num_samples_per_call 20 
#                           --reference /path/to/reference/genome 
#                           --tiebrush /path/to/tiebrush/executable 
#                           --tiewrap /path/to/tiewrap/executable 
#                           --sample_counts
# - input csv should be in the format of:
#       Tissue,/path/to/sample
#       Tissue,/path/to/sample
#       ...
#===================================================================

import os
import sys
import argparse
import subprocess

def run_tissue_tiebrush(args,t2s):
    outdir = args.outdir.rstrip("/")+"/"
    tb_fname = outdir+"tb.txt"
    with open(tb_fname,"w+") as tbFP:
        for tissue,paths in t2s.items():
            tissue_dir=outdir+tissue+"/"
            if not os.path.exists(tissue_dir):
                os.makedirs(tissue_dir)
            
            lst_fname = outdir+tissue+".lst"
            with open(lst_fname,"w+") as lstFP:
                lstFP.write("\n".join(paths))
            
            tbFP_output = args.tiewrap+" --batch-size "+str(args.num_samples_per_call)+" --output "+tissue_dir+tissue+".tb.bam"+" "+" "+lst_fname+"\n"
            if (args.sample_counts):
                tbFP_output = args.tiewrap+" --sample-counts --batch-size "+str(args.num_samples_per_call)+" --output "+tissue_dir+tissue+".tb.bam"+" "+" "+lst_fname+"\n"

            tbFP.write(tbFP_output)
            
            
    # now can run the file in parallel
    parallel_cmd = "parallel -j "+str(args.threads)+" < "+tb_fname
    subprocess.call(parallel_cmd,shell=True)
    
    print("indexing tissue merges")
    for tissue,paths in t2s.items():
        tissue_dir=outdir+tissue+"/"
        idx_cmd = ["samtools","index",tissue_dir+tissue+".tb.bam"]
        subprocess.call(idx_cmd)    

    print("merging all tmps")

    # lastly run the final tiebrush to merge them all together
    tiewrap_cmd = [args.tiewrap,
                   "--batch-size",str(args.num_samples_per_call),
                   "--output",outdir+"all.tb.bam"]
    if (args.sample_counts):
        tiewrap_cmd = [args.tiewrap,
                       "--sample-counts",
                       "--batch-size",str(args.num_samples_per_call),
                       "--output",outdir+"all.tb.bam"]
    for tissue,paths in t2s.items():
        tissue_dir=outdir+tissue+"/"
        tiewrap_cmd.append(tissue_dir+tissue+".tb.bam")
    print(" ".join(tiewrap_cmd))
    subprocess.call(tiewrap_cmd)

    idx_cmd = ["samtools","index",outdir+"all.tb.bam"]

def run_tissue_tiecov_default(args,t2s):
    outdir = args.outdir.rstrip("/")+"/"
    with open(outdir+"tiecov_default.parallel","w+") as outFP:
        for tissue,paths in t2s.items():
            print("run_tissue_tiecov_default: "+tissue)
            tissue_dir = outdir+tissue+"/"
            outFP.write(args.tiecov_default+" -s "+tissue_dir+tissue+".def.sample -j "+tissue_dir+tissue+".def.junctions -c "+tissue_dir+tissue+".def.coverage "+tissue_dir+tissue+".tb.bam\n")

    # now can run the file in parallel
    parallel_cmd = "parallel -j "+str(args.threads)+" < "+outdir+"tiecov_default.parallel"
    subprocess.call(parallel_cmd,shell=True)

def run_convert_to_bigwig(args,t2s):
    outdir = args.outdir.rstrip("/")+"/"
    with open(outdir+"bigwig.parallel","w+") as outFP:
        for tissue,paths in t2s.items():
            print("run_convert_to_bigwig: "+tissue)
            tissue_dir = outdir+tissue+"/"
            outFP.write("bedtools sort -i "+tissue_dir+tissue+".def.coverage.bedgraph > "+tissue_dir+tissue+".def.coverage.sorted.bedgraph && bedGraphToBigWig "+tissue_dir+tissue+".def.coverage.sorted.bedgraph ~/genomes/human/hg38/hg38_p12_ucsc.no_alts.no_fixs.fa.fai "+tissue_dir+tissue+".def.coverage.bigwig\n")

    # now can run the file in parallel
    parallel_cmd = "parallel -j "+str(args.threads)+" < "+outdir+"bigwig.parallel"
    
    subprocess.call(parallel_cmd,shell=True)

def run_tiebrush_tiecov_bigwig_all(args,t2s):
    print("run_tiebrush_tiecov_bigwig_all")
    outdir = args.outdir.rstrip("/")+"/"

    tiecov_cmd = [args.tiecov_default,
                "-s",outdir+"all.def.sample",
                "-j",outdir+"all.def.junctions",
                "-c",outdir+"all.def.coverage",outdir+"all.tb.bam"]
    print(" ".join(tiecov_cmd))
    subprocess.call(tiecov_cmd)

    sort_cmd =  "bedtools sort -i "+outdir+"all.def.coverage.bedgraph > "+outdir+"all.def.coverage.sorted.bedgraph"
    print(sort_cmd)
    subprocess.call(sort_cmd,shell=True)

    bw_cmd =  ["bedGraphToBigWig",
                 outdir+"all.def.coverage.sorted.bedgraph",args.reference,outdir+"all.def.coverage.bigwig"]
    print(" ".join(bw_cmd))
    subprocess.call(bw_cmd)

def run_step1(args):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    assert os.path.exists(args.input),"input file does not exist"

    # first need to form a dictionary of tissues to samples
    tissue2samples = dict()
    with open(args.input,"r") as inFP:
        for line in inFP.readlines():
            line = line.strip()
            tissue,cram_fp = line.split(",")
            tissue2samples.setdefault(tissue,[]).append(cram_fp)

    run_tissue_tiebrush(args,tissue2samples)
    run_tissue_tiecov_default(args,tissue2samples)
    run_convert_to_bigwig(args,tissue2samples)
    run_tiebrush_tiecov_bigwig_all(args,tissue2samples)

def main(args):
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('--input',
                        required=True,
                        type=str,
                        help="Input file in CSV format where column #1 is the path to a cram file and column #2 is the tissue name")
    parser.add_argument('--outdir',
                        required=True,
                        type=str,
                        help="Output directory in which all output and temporary data will be stored")
    parser.add_argument("--threads",
                        required=False,
                        type=int,
                        default=1,
                        help="number of threads to be used by GNU parallel")
    parser.add_argument("--num_samples_per_call",
                        required=False,
                        type=int,
                        default=20,
                        help="number of samples to process with tiebrush within a single batch")
    parser.add_argument("--tiebrush",
                        required=False,
                        type=str,
                        default="tiebrush",
                        help="path to the tiebrush executable")
    parser.add_argument("--tiewrap",
                        required=False,
                        type=str,
                        default="tiewrap.py",
                        help="path to the tiewrap executable")
    parser.add_argument("--tiecov_default",
                        required=False,
                        type=str,
                        default="tiecov",
                        help="path to the standard tiecov executable")
    parser.add_argument("--reference",
                        required=True,
                        type=str,
                        help="path to the reference genome")
    parser.add_argument("--sample_counts",
                        required=False,
                        default=False,
                        action='store_true',
                        help="toggle sample count tracking")
    

    parser.set_defaults(func=run_step1)
    args=parser.parse_args()
    args.func(args)


if __name__=="__main__":
    main(sys.argv[1:])
