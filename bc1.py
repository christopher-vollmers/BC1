### written by Christopher Vollmers ###

import argparse
import mappy
import os
import sys
import gc
from time import localtime, strftime


PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
sys.path.append(os.path.abspath(PATH))

BD1Path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
abpoa=BD1Path+'/abPOA-v1.4.1/bin/abpoa'
racon=BD1Path+'/racon/build/bin/racon'
minimap2=BD1Path+'/minimap2/minimap2'

VERSION = "v0.95 - The readings are off the charts. Over Q50. Even Master Yoda doesn't have a Midi-chlorian count that high"

from utils import extracting_UMIs,extracting_UMIs_STARsolo,create_labeled_subreads,collect_indexes,find_fuzzy_matches,write_new_indexes,read_fasta,findUMISequence,read_fasta_root
from utils import delegating_consensus_generation, process_batch, determine_consensus,read_cons,write_chimeras,split_reads


parser=argparse.ArgumentParser()
parser.add_argument('-i','--input_fasta',type=str, help='demultiplexed,trimmed, and reoriented R2C2 consensus reads as output by C3POa_postprocessing.py')
parser.add_argument('-o','--output_file_root',type=str, help='root file name. Intermediate and final output files will be start with this')
parser.add_argument('-p','--output_path',type=str, help='output directory. Intermediate and final output files will be in this directory')
parser.add_argument('-r','--resume',action='store_true',default=False, help='if set, bc1 will try to resume from previous output')
parser.add_argument('-s','--subread_files',type=str, help='R2C2 subreads as output by C3POa. Can be multiple comma separated files')
parser.add_argument('-m','--medaka',action='store_true',default=False, help='if set, medaka will be run. required for highest accuracy but definitely slows things down. assumes data was produced with R10.4 pores at 400bp/s speed')
parser.add_argument('-f','--fuzzy',action='store_true',default=False, help='if set, reads with a single mismatch between their UMIs will be combined. By default, only reads with identical UMIs will be combined.')
#parser.add_argument('-S','--STARsolo_sam_file',type=str, help='experimental setting.Probably best left alone for now')
parser.add_argument('-b','--subsample',type=int,default=1000,help='subreads will be subsampled to this number.' )
parser.add_argument('-t', '--threads', type=int, help='defines the number of threads the multiprocessing will use')
parser.add_argument('-u','--UMIpatterns',type=str, help="""UMI patterns separated by commas. \nAn example would be '5.0:3.GACAG.NNNNNNNNNNNNNN.,3.0:3.CAC.NNNNNN.TTTT' that would indicate two UMIs.
                                                         The first at the 5prime end of the read starting somewhere in the first 3 bases of the read and is flanked on the left with 'GACAG' and is 14nt long.
                                                         The second at the 3prime end of the read starting somewhere in the first 3 bases of the (reverse complemented) read and is flanked on the left with 'CAC' and on the right with 'TTTT' and is 6nt long. IUPAC wild card bares can be used at any position""")



args=parser.parse_args()
output_file_root=args.output_file_root
output_path=os.path.abspath(args.output_path)
input_reads = os.path.abspath(args.input_fasta)
subread_file_prelim=args.subread_files
subreads=[]
if subread_file_prelim:
    for subread_files in subread_file_prelim.split(','):
        subreads.append(os.path.abspath(subread_files))
else:
    print('no subreads provided. Will use regular input reads for error correction') 
    subreads.append(input_reads)
subread_file=(',').join(subreads)

threads=args.threads
fuzzy=args.fuzzy
medaka=args.medaka
subsample=args.subsample
resume=args.resume




UMIpatterns=args.UMIpatterns
combined_umi_length=0
for UMIpattern in UMIpatterns.split(','):
    combined_umi_length+=len(UMIpattern.split('.')[3])
combined_umi_length+=(len(UMIpatterns.split(','))-1)

def main():

    if not os.path.isdir(output_path):
        os.system(f'mkdir {output_path}')
    if not os.path.isdir(f'{output_path}/tmp'):
        os.system(f'mkdir {output_path}/tmp')

 #   extractedUMIs=False
    labeledSubreads=False
    matchedFuzzy=False
    sortedSubreads=False
    splitSubreads=False
    processed=set()


    print('\nBC1 "'+ VERSION +'" was run on '\
                  + strftime("%Y-%m-%d %H:%M:%S", localtime())+'\n'\
                  + 'with the following parameters\n'\
                  + str(args).replace('Namespace(','').replace(')','') + '\n')


    logfile=f'{output_path}/{output_file_root}.log'

    if not resume:
        logfile_object=open(logfile,'w')
    else:
        print('resuming previous run')
        for line in open(logfile):
            if line.strip()==f'labeled subreads\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}':
                labeledSubreads=True
            if line.strip()==f'matched fuzzy\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}':
                matchedFuzzy=True
            if line.strip()==f'sorted subreads\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}':
                sortedSubreads=True
            if line.strip()==f'split subreads\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}':
                splitSubreads=True
            if f'processed\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}' in line.strip():
                processed.add(line.strip().split('\t')[0])

        print(f'labeling subreads completed:{labeledSubreads}\nfuzzy matching of UMIs completed (optional):{matchedFuzzy}\nsorting subreads completed:{sortedSubreads}\nsplitting subreads completed:{splitSubreads}\nfiles already processed:{processed}\n')

        logfile_object=open(logfile,'a')

    logfile_object.write('\nBC1 "'+ VERSION +'" was run on '\
                  + strftime("%Y-%m-%d %H:%M:%S", localtime())+'\n'\
                  + 'with the following parameters\n'\
                  + str(args).replace('Namespace(','').replace(')','') + '\n')
    logfile_object.flush()
    if not labeledSubreads:
        print('extracting UMIs from', input_reads)
        UMIdict = extracting_UMIs(input_reads,f'{output_path}/tmp/{output_file_root}',UMIpatterns)
        logfile_object.write(f'extracted UMIs\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}\n')
        logfile_object.flush()
        print('labeling subreads in', subread_file, '(this can take a while)')
        create_labeled_subreads(UMIdict,subread_file,f'{output_path}/tmp/{output_file_root}.subreads')
        logfile_object.write(f'labeled subreads\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}\n')
        logfile_object.flush()
    toSort=f'{output_path}/tmp/{output_file_root}.subreads'

    if fuzzy:
        if not matchedFuzzy:
            all_indexes,first_indexes,second_indexes = collect_indexes(f'{output_path}/tmp/{output_file_root}.subreads')
            print('combine similar UMIs and their Indexes')
            equivalent_indexes,chimeras = find_fuzzy_matches(all_indexes,first_indexes,second_indexes)
            print('writing and sorting subreads')
            write_new_indexes(f'{output_path}/tmp/{output_file_root}.subreads',equivalent_indexes)
            logfile_object.write(f'matched fuzzy\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}\n')
            logfile_object.flush()
        toSort=f'{output_path}/tmp/{output_file_root}.subreads.fuzzy'

    if not sortedSubreads:
        print('sorting subreads')
        os.system(f'sort --parallel {threads} -k1,1n -k2,2n -T {output_path} {toSort} > {output_path}/tmp/{output_file_root}.subreads.sorted')
        logfile_object.write(f'sorted subreads\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}\n')
        logfile_object.flush()
    if not splitSubreads:
         split_reads(f'{output_path}/tmp/{output_file_root}.subreads.sorted',input_reads,f'{output_path}/tmp/{output_file_root}.UMIs',1000000)
         logfile_object.write(f'split subreads\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}\n')
         logfile_object.flush()

    split_files=[]
    for read_file in os.listdir(f'{output_path}/tmp/'):
        absFile=f'{output_path}/tmp/{read_file}'
        if os.path.isfile(absFile):
            if f'{output_file_root}.subreads.sorted.split.' in read_file and 'merged' not in read_file:
                split_files.append(absFile)
                if absFile not in processed:
                    print(read_file)
                    if not medaka:
                        os.system(f'python3 {BD1Path}/createConsensusSequences.py -t {output_path}/tmp/ -i {absFile} -u {combined_umi_length} -n {threads} -b {subsample}')
                    else:
                        os.system(f'python3 {BD1Path}/createConsensusSequences.py -t {output_path}/tmp/ -i {absFile} -u {combined_umi_length} -n {threads} -b {subsample} -m')
                    logfile_object.write(f'{absFile}\tprocessed\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}\n')
                    logfile_object.flush()
                else:
                    print(f'skipping {absFile}. Already processed')

    print('merging output files')
    outMerged=open(f'{output_path}/{output_file_root}.merged.fasta','w')
    for absFile in split_files:
        mergedFile=f'{absFile}.merged.fasta'
        for name,seq,q in mappy.fastx_read(mergedFile):
            outMerged.write(f'>{name}\n{seq}\n')
    outMerged.close()

    logfile_object.write(f'merged output\t{output_path}\t{output_file_root}\t{input_reads}\t{subread_file}\t{UMIpatterns}\t{fuzzy}\t{medaka}\t{subsample}\n')
    logfile_object.close()
main()
