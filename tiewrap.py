#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import multiprocessing

def build_base_cmd(args):
    # check if tiebrush found
    tiebrush_path = "/".join(os.path.realpath(__file__).split("/")[:-1])+"/"+"tiebrush"
    assert os.path.exists(tiebrush_path),"tiebrush executable not found. please make sure it is located in the same directory as tiewrap.py"

    tb_cmd = [tiebrush_path]
    if args.cigar:
        tb_cmd.append("-C")
    if args.clip:
        tb_cmd.append("-P")
    if args.exon:
        tb_cmd.append("-E")
    if args.keep_unmap:
        tb_cmd.append("-M")
    if args.keep_supp:
        tb_cmd.append("-S")
    if args.max_nh:
        tb_cmd.append("-N")
        tb_cmd.append(str(args.max_nh))
    if args.max_map_qual:
        tb_cmd.append("-Q")
        tb_cmd.append(str(args.max_map_qual))
    if args.flags:
        tb_cmd.append("-F")
        tb_cmd.append(str(args.flags))

    return " ".join(tb_cmd)

def create_batches(fnames,tb_cmd,args,rep):
    batches = []
    batch_fnames = []
    for i in range(0,len(fnames),args.batch_size):
        out_fname = args.output+".b"+str(rep)+"."+str(i)+".bam"
        batches.append(tb_cmd+" -o "+out_fname+" "+" ".join(fnames[i:i + args.batch_size]))
        batch_fnames.append(out_fname)

    return batches,batch_fnames

def run(cmd):
    cmd = cmd.split(" ")
    subprocess.call(cmd)

def start(args):
    tb_cmd = build_base_cmd(args)

    # first determine whether inputs are listed in a file
    # or they are listed one-by-one at the end
    assert len(args.rest)>0,"no input files provided"

    # check if all sam/bam/cram files
    sam=True
    for fname in args.rest:
        if fname.split(".")[-1] in ["sam","bam","cram"]:
            if not sam:
                print("check your inputs")
                exit(1)
            sam=True
        else:
            sam=False

    if not sam and len(args.rest)>1:
        print("check your inputs")
        exit(1)

    fnames = []
    if sam:
        fnames = args.rest
    else:
        assert len(args.rest)==1,"check your inputs"
        with open(args.rest[0],"r") as inFP:
            for line in inFP.readlines():
                line = line.strip()
                assert line.split(".")[-1] in ["sam","bam","cram"],"unrecognized file format: "+line
                fnames.append(line)

    # next check that all files are valid (exist)
    for fname in fnames:
        assert os.path.exists(fname),"file does not exist: "+fname

    # next batch them if user requests batching or multithreaded
    round=0
    batches,batch_fnames = create_batches(fnames,tb_cmd,args,round)
    round+=1

    tmp_fnames = [] # running list of files to be removed

    if len(batches)>1:
        # run batches
        while(len(batches)>1):
            with multiprocessing.Pool(processes=args.threads) as pool:
                pool.map(run,batches)

            for i in range(len(tmp_fnames)):
                assert os.path.exists(tmp_fnames[i]),"tmp file does not exist and cannot be removed: "+tmp_fnames[i]
                os.remove(tmp_fnames[i])

            tmp_fnames = []
            tmp_fnames.extend(batch_fnames)

            # repeat batching until done
            old_batch_fnames = batch_fnames
            batches,batch_fnames = create_batches(batch_fnames,tb_cmd,args,round)
            round+=1
        batch_fnames=old_batch_fnames # since the last set did not get executed within the while loop

    else: # a single batch - run final
        batch_fnames = fnames

    # run final merge
    batches = tb_cmd+" -o "+args.output+" "+" ".join(batch_fnames)
    run(batches)

    # cleanup any tmp data that might still exist
    for i in range(len(tmp_fnames)):
        assert os.path.exists(tmp_fnames[i]),"tmp file does not exist and cannot be removed: "+tmp_fnames[i]
        os.remove(tmp_fnames[i])

    return

def tiebrush(argv):
    parser = argparse.ArgumentParser(description='''Help Page''')

    parser.add_argument('-o',
                        '--output',
                        required=True,
                        type=str,
                        help="output file")
    parser.add_argument('-C',
                        '--cigar',
                        required=False,
                        action='store_true',
                        help="merge if only CIGAR string is the same")
    parser.add_argument('-P',
                        '--clip',
                        required=False,
                        action='store_true',
                        help="merge if clipped CIGAR string is the same")
    parser.add_argument('-E',
                        '--exon',
                        required=False,
                        action='store_true',
                        help="merge if exon boundaries are the same")
    parser.add_argument('-S',
                        '--keep-supp',
                        required=False,
                        action='store_true',
                        help="keep supplementary alignments")
    parser.add_argument('-M',
                        '--keep-unmap',
                        required=False,
                        action='store_true',
                        help="keep unmapped reads")
    parser.add_argument('-N',
                        '--max-nh',
                        required=False,
                        type=int,
                        help="maximum NH score of the reads to retain")
    parser.add_argument('-Q',
                        '--max-map-qual',
                        required=False,
                        type=int,
                        help="maximum NH score of the reads to retain")
    parser.add_argument('-F',
                        '--flags',
                        required=False,
                        type=int,
                        help="bits in SAM flag to use in read comparison")
    parser.add_argument('-t',
                        '--threads',
                        required=False,
                        default=1,
                        type=int,
                        help="number of threads to use")
    parser.add_argument('-b',
                        '--batch-size',
                        required=False,
                        default=100,
                        type=int,
                        help="Number of input files to process in a batch on each thread")
    parser.add_argument('rest', nargs=argparse.REMAINDER)

    parser.set_defaults(func=start)

    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
    tiebrush(sys.argv[1:])