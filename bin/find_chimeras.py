import sys
import mappy as mm
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--inFile', '-i', type=str, action='store', help='file to be filtered')
parser.add_argument('--out_root', '-o', type=str, action='store', help='output file root')
parser.add_argument('--reference', '-r', type=str, action='store', help='fasta file containing reference sequences')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(0)
args = parser.parse_args()
inFile = args.inFile
out_root = args.out_root
reference=args.reference


mm_align = mm.Aligner(reference, preset='sr',best_n=50)
if not mm_align: raise Exception("ERROR: failed to load/build index")

all=set()
for name,seq,q in mm.fastx_read(reference):
    all.add(name)

out1=open(f'{out_root}.chimera.fasta','w')
out2=open(f'{out_root}.nochimera.fasta','w')

for name,seq,q in mm.fastx_read(inFile):
#  show=False
#  if 'medaka' in name and int(name.split('_')[6])>=4:
  show=True
  if show:
    print(name)
    switch=0
    potential=all
    for pos in range(0,len(seq),10):
        partSeq=seq[pos:pos+50]
        match=''
        possibles=[]
        possibleSet=set()
        for hit in mm_align.map(partSeq):
            possibles.append((hit.NM,hit.ctg,pos,hit.r_st))
        possibles=sorted(possibles)
        if possibles:
            lowest=possibles[0][0]
            if lowest==0:
                for possible in possibles:
                    if possible[0]==lowest:
                        possibleSet.add(possible[1])
                previous=potential
                potential=potential&possibleSet
                if len(potential)==0:
                    switch+=1
                    potential=possibleSet


    if switch>0:
        out1.write(f'>{name}\n{seq}\n')
    else:
        out2.write(f'>{name}\n{seq}\n')
