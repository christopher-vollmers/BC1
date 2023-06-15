import sys
import mappy
import editdistance

primers=sys.argv[1]
consensus_reads=sys.argv[2]
out=open(sys.argv[3],'w')

leftUMIlength=58
rightUMIlength=58

IUPACdict={}
IUPACdict['A']=['A']
IUPACdict['T']=['T']
IUPACdict['G']=['G']
IUPACdict['C']=['C']
IUPACdict['R']=['A','G']
IUPACdict['Y']=['C','T']
IUPACdict['S']=['G','C']
IUPACdict['W']=['A','T']
IUPACdict['K']=['G','T']
IUPACdict['M']=['A','C']
IUPACdict['B']=['C','G','T']
IUPACdict['D']=['A','G','T']
IUPACdict['H']=['A','C','T']
IUPACdict['V']=['A','C','G']
IUPACdict['N']=['A','T','G','C']


primerDict={}
for name,seq,q in mappy.fastx_read(primers):
    primerDict[name]=[]
    seq=seq.upper().strip()
    seqs=['']
    for base in seq:
        newSeqs=[]
        print(IUPACdict[base])
        for unambigousBase in IUPACdict[base]:
            print(unambigousBase)
            for unambigousSeq in seqs:
                newSeqs.append(unambigousSeq+unambigousBase)
        seqs=newSeqs
        print(seqs)
    for seq in seqs:
        primerDict[name].append(seq)

print(primerDict)


for name,seq,q in mappy.fastx_read(consensus_reads):
    forwardDists=[]
    reverseDists=[]
    forwardSeq=seq
    reverseSeq=mappy.revcomp(seq)
    for primer,seqs in primerDict.items():
        for seq2 in seqs:
            for wiggle in range(-6,6,1):
                subseq1=forwardSeq[leftUMIlength+wiggle:leftUMIlength+wiggle+len(seq2)]
                subseq2=reverseSeq[rightUMIlength+wiggle:rightUMIlength+wiggle+len(seq2)]
                dist1=editdistance.eval(subseq1,seq2)
                dist2=editdistance.eval(subseq2,seq2)
                forwardDists.append((dist1,seq2,subseq1,primer,leftUMIlength+wiggle+len(seq2)))
                reverseDists.append((dist2,seq2,subseq2,primer,rightUMIlength+wiggle+len(seq2)))
    bestForward=sorted(forwardDists)[0]
    bestReverse=sorted(reverseDists)[0]

    if bestForward[0]==0 and bestReverse[0]==0:
       out.write('>%s\n%s\n' % (name+'|'+bestForward[3]+'|'+bestReverse[3],seq[bestForward[4]:(bestReverse[4])*-1]))
