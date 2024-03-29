### written by Christopher Vollmers ###

import argparse
import mappy
import os
import sys
import gc

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
sys.path.append(os.path.abspath(PATH))

BD1Path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
abpoa=BD1Path+'/abPOA-v1.4.1/bin/abpoa'
racon=BD1Path+'/racon/build/bin/racon'
minimap2=BD1Path+'/minimap2/minimap2'

from utils import process_batch


parser=argparse.ArgumentParser()
parser.add_argument('-t','--temp_folder',type=str, help='temp folder')
parser.add_argument('-i','--input_file',type=str, help='File generated by the BC1 main script')
parser.add_argument('-m','--medaka',action='store_true',default=False, help='if set, medaka will be run. required for highest accuracy but definitely slows things down. assumes data was produced with R10.4 pores at 400bp/s speed')
parser.add_argument('-b','--subsample',type=int,default=1000,help='subreads will be subsampled to this number.' )
parser.add_argument('-n', '--threads', type=int, help='defines the number of threads the multiprocessing will use')
parser.add_argument('-u', '--combined_umi_length', type=int, help='combined length (including spacers) of all the UMIs in the reads')


args=parser.parse_args()
input_file=args.input_file
combined_UMI_length=args.combined_umi_length
threads=args.threads
medaka=args.medaka
subsample=args.subsample
temp_folder=args.temp_folder


def main():

    singles = 0
    processed = 0

    out=open(f'{input_file}.merged.fasta','w')
    subreads={}
    reads={}
    for line in open(input_file):
        a=line.strip().split('\t')
        type1,UMInumber,UMI,name,seq,qual = a[0],a[1],a[2],a[3],a[4],a[5]
        if UMInumber not in subreads:
            subreads[UMInumber]=set()
            reads[UMInumber]=set()
        if type1=='cons':
            reads[UMInumber].add((name,seq,qual,UMI))
        elif type1=='sub':
            subreads[UMInumber].add((name,seq,qual))

    processed,singles=process_batch(subreads,reads,processed,singles,subsample,medaka,threads,out,combined_UMI_length,abpoa,racon,minimap2,temp_folder)
    out.close()
    print(f'{processed} reads merged from multi-read UMIs, {singles} reads not merged because incomplete or single-read UMI', ' '*60)


main()
