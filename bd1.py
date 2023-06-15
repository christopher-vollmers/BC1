### written by Christopher Vollmers ###

import argparse
import mappy
import os
import sys
import gc

PATH = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/bin/'
sys.path.append(os.path.abspath(PATH))

from utils import extracting_UMIs,extracting_UMIs_STARsolo,create_labeled_subreads,collect_indexes,find_fuzzy_matches,write_new_indexes,read_fasta,findUMISequence,read_fasta_root
from utils import delegating_consensus_generation, process_batch, determine_consensus,read_cons,write_chimeras


parser=argparse.ArgumentParser()
parser.add_argument('-i','--input_fasta',type=str)
parser.add_argument('-o','--output_file_root',type=str)
parser.add_argument('-s','--subread_files',type=str)
parser.add_argument('-m','--medaka',action='store_true',default=False)
parser.add_argument('-f','--fuzzy',action='store_true',default=False)
parser.add_argument('-S','--STARsolo_sam_file',type=str)
parser.add_argument('-b','--subsample',type=int,default=200)
parser.add_argument('-t', '--threads', type=int, help='defines the number of threads the multiprocessing will use')
parser.add_argument('-u','--UMIpatterns',type=str, help="""UMI patterns separated by commas. \nAn example would be '5.0:3.GACAG.NNNNNNNNNNNNNN.,3.0:3.CAC.NNNNNN.TTTT' that would indicate two UMIs.
                                                         The first at the 5prime end of the read starting somewhere in the first 3 bases of the read and is flanked on the left with 'GACAG' and is 14nt long.
                                                         The second at the 3prime end of the read starting somewhere in the first 3 bases of the (reverse complemented) read and is flanked on the left with 'CAC' and on the right with 'TTTT' and is 6nt long""")



args=parser.parse_args()
output_file_root=args.output_file_root
input_reads = args.input_fasta
subread_file=args.subread_files
UMIpatterns=args.UMIpatterns
threads=args.threads
fuzzy=args.fuzzy
medaka=args.medaka
subsample=args.subsample
STARsolo=args.STARsolo_sam_file


def main():
    if STARsolo:
        print('extracting UMIs from', STARsolo)
        UMIdict,combined_UMI_length = extracting_UMIs_STARsolo(STARsolo,output_file_root,UMIpatterns)
    else:
        print('extracting UMIs from', input_reads)
        UMIdict,combined_UMI_length = extracting_UMIs(input_reads,output_file_root,UMIpatterns)
    combined_UMI_length = 36
    print('combined UMI pattern length',combined_UMI_length)
    print('labeling and sorting subreads in', subread_file, '(this can take a while)')
    create_labeled_subreads(UMIdict,subread_file,output_file_root+'.subreads')
    print('UMI collecting indexes')
    toSort=f'{output_file_root}.subreads'
    if fuzzy:
        all_indexes,first_indexes,second_indexes = collect_indexes(output_file_root+'.subreads')
        print('combine similar UMIs and their Indexes')
        equivalent_indexes,chimeras = find_fuzzy_matches(all_indexes,first_indexes,second_indexes)
        print('writing and sorting subreads')
        write_new_indexes(output_file_root+'.subreads',equivalent_indexes)
        write_chimeras(output_file_root+'.chimeras',chimeras)
        toSort=f'{output_file_root}.subreads.fuzzy'
    os.system(f'sort --parallel {threads} -k1,1n -k2,2n -T ./ {toSort} > {output_file_root}.subreads.sorted')
    out=open(output_file_root+'.merged.fasta','w')
    delegating_consensus_generation(input_reads,output_file_root+'.UMIs',output_file_root+'.subreads.sorted',subsample,medaka,combined_UMI_length,threads,1000000,out)


main()
