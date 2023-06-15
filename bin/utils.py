import sys
import mappy
import numpy as np
import editdistance
import multiprocessing as mp
import pyabpoa as poa
import gc
import os

poa_aligner = poa.msa_aligner(match=5)


def create_labeled_subreads(UMIdict,subread_files,out2):
    outSub=open(out2,'w')
    for subread_file in subread_files.split(','):
        for name,seq,q in mappy.fastx_read(subread_file):
            root=name.split('_')[0]
            if root in UMIdict:
                UMI,number=UMIdict[root]
                outSub.write('%s\t%s\t%s\t%s\t%s\n' %(str(number),UMI,name,seq,q))
    outSub.close()


def findUMISequence(sequence,pattern):

    IUPACdict={}
    IUPACdict['A']=set(['A'])
    IUPACdict['T']=set(['T'])
    IUPACdict['G']=set(['G'])
    IUPACdict['C']=set(['C'])
    IUPACdict['R']=set(['A','G'])
    IUPACdict['Y']=set(['C','T'])
    IUPACdict['S']=set(['G','C'])
    IUPACdict['W']=set(['A','T'])
    IUPACdict['K']=set(['G','T'])
    IUPACdict['M']=set(['A','C'])
    IUPACdict['B']=set(['C','G','T'])
    IUPACdict['D']=set(['A','G','T'])
    IUPACdict['H']=set(['A','C','T'])
    IUPACdict['V']=set(['A','C','G'])
    IUPACdict['N']=set(['A','T','G','C'])
    UMI=''
    valid=True
    direction,range1,left,variable,right = pattern.split('.')
    if direction not in ['5','3']:
        print('invalid pattern, direction has to be 5 or 3')
        valid=False
        reason='invalid direction'
    if direction=='3':
        sequence=mappy.revcomp(sequence)

    start,end = int(range1.split(':')[0]),int(range1.split(':')[1])
    left_variable_start=''
    right_start=''
    if len(sequence)<100:
        valid=False
        reason='sequence shorter than 100nt'

    if valid:
        UMIpattern=left+variable+right
        valid=False
        reason='no UMI pattern match'
        for pos in range(start,end,1):
            matches=0
            match=sequence[pos:pos+len(UMIpattern)].upper()
            for index in range(0,len(UMIpattern),1):
                v_base=UMIpattern[index]
                s_base=match[index]
                if s_base in IUPACdict[v_base]:
                    matches+=1
            if len(UMIpattern)==matches:
                if right:
                    UMI=match[len(left):-len(right)]
                else:
                    UMI=match[len(left):]
                valid=True
                break
    if valid:
        return UMI,'UMI found that matched pattern '+pattern
    else:
        return '', reason

def read_fasta_root(infile):

    reads={}
    for name,seq,qual in mappy.fastx_read(infile):
        name=name.split()[0]
        name_root=name.split('_')[0]
        reads[name_root]=(name,seq)

    return reads


def extracting_UMIs(input_reads,output_file_root,UMIpatterns):

    reads=read_fasta_root(input_reads)
    out1=open(output_file_root+'.UMIs','w')
    UMIdict={}
    UMInumbers={}
    counter=0
    reasonDict={}
    for name_root,(name,sequence) in reads.items():
        UMIs=[]
        for UMIpattern in UMIpatterns.split(','):
            UMI,reason=findUMISequence(sequence,UMIpattern)
            if reason not in reasonDict:
                reasonDict[reason]=0
            reasonDict[reason]+=1
            UMIs.append(UMI)
        UMI=('_').join(UMIs)
        if UMI not in UMInumbers:
            counter+=1
            UMInumbers[UMI]=counter

        out1.write('%s\t%s\n' %(name,UMI))
        root=name.split('_')[0]
        UMIdict[root]=(UMI,UMInumbers[UMI])

    print(reasonDict.items())
    combined_umi_length=0
    for UMIpattern in UMIpatterns.split(','):
        combined_umi_length+=len(UMIpattern.split('.')[3])
    combined_umi_length+=(len(UMIpatterns.split(','))-1)
    return UMIdict,combined_umi_length


def extracting_UMIs_STARsolo(input_sam,output_file_root,UMIpatterns):

    out1=open(output_file_root+'.UMIs','w')
    UMIdict={}
    UMInumbers={}
    counter=0
    reasonDict={}
    lengths=[]
    for line in open(input_sam):
        UMIs=[]
        a=line.strip().split('\t')
        name=a[0]
        for UMIpattern in UMIpatterns.split(','):
            tags=a[8:]
            UMI=''
            for tag in tags:
                if tag.startswith(UMIpattern):
                    UMI=tag.split(UMIpattern)[1].replace('_','')
                    break
            if not UMI:
                reason='no '+UMIpattern+' found'
            else:
                reason=UMIpattern+' found'
            if reason not in reasonDict:
                reasonDict[reason]=0
            reasonDict[reason]+=1
            UMIs.append(UMI)

        UMI=('_').join(UMIs)
        if UMI not in UMInumbers:
            counter+=1
            UMInumbers[UMI]=counter
            lengths.append(len(UMI))
        out1.write('%s\t%s\n' %(name,UMI))
        root=name.split('_')[0]
        UMIdict[root]=(UMI,UMInumbers[UMI])

    print(reasonDict.items())
    combined_umi_length=int(np.median(lengths))
    return UMIdict,combined_umi_length

def collect_indexes(inFile):
    counter=0
    first_indexes={}
    second_indexes={}
    all_indexes=set()
    for line in open(inFile):
        counter+=1
        a=line.strip().split('\t')
        index_number=a[0]
        index=a[1]
        if '_' in index:
            index_sequence1,index_sequence2=index.split('_')[0],index.split('_')[1]
        else:
            length=len(index)
            index_sequence1,index_sequence2=index[:int(length/2)],index[int(length/2):]

        if index_sequence1 not in first_indexes:
            first_indexes[index_sequence1]={}
        first_indexes[index_sequence1][index_sequence2]=index_number

        if index_sequence2 not in second_indexes:
            second_indexes[index_sequence2]={}
        second_indexes[index_sequence2][index_sequence1]=index_number
        all_indexes.add((index_sequence1,index_sequence2,index_number))

    return all_indexes, first_indexes, second_indexes


def find_fuzzy_matches(all_indexes,first_indexes,second_indexes):
    equivalent_indexes={}
    chimeras=set()
    sortedIndexes=sorted(list(all_indexes),key=lambda x:int(x[2]))
    lastIndex=sortedIndexes[-1]
    for index in sortedIndexes:
        print('working on index number', index, 'of',lastIndex,end='\r')
        index_number=index[2]
        index_sequence1=index[0]
        index_sequence2=index[1]
        if index_number in equivalent_indexes:
            new_number=equivalent_indexes[index_number]
        else:
            new_number=index_number

        for new_index in first_indexes[index_sequence1]:
            compare_index_number=first_indexes[index_sequence1][new_index]
            if compare_index_number!=index_number:
                if len(index_sequence2)==len(new_index):
                    dis1=editdistance.eval(index_sequence2,new_index)
                    if dis1<=1:
                        equivalent_indexes[compare_index_number]=new_number
                    else:
                        chimeras.add(compare_index_number)
                        chimeras.add(new_number)

        for new_index in second_indexes[index_sequence2]:
            compare_index_number=second_indexes[index_sequence2][new_index]
            if compare_index_number!=index_number:
                if len(index_sequence1)==len(new_index):
                    dis1=editdistance.eval(index_sequence1,new_index)
                    if dis1<=1:
                        equivalent_indexes[compare_index_number]=new_number
                    else:
                        chimeras.add(compare_index_number)
                        chimeras.add(new_number)

    print('merged indexes',len(equivalent_indexes),'likely chimeras',len(chimeras))
    return equivalent_indexes,chimeras

def write_new_indexes(inFile,equivalent_indexes):
    out=open(inFile+'.fuzzy','w')
    for line in open(inFile):
        a=line.strip().split('\t')
        index_number=a[0]
        if index_number in equivalent_indexes:
            a[0]=equivalent_indexes[index_number]+'M'
        new_line=('\t').join(a)+'\n'
        out.write(new_line)
    out.close()

def write_chimeras(outFile,chimeras):
    out=open(outFile,'w')
    for chimera in chimeras:
        out.write(chimera+'\n')
    out.close()


def read_fasta(infile):
    reads = {}
    for name,seq,qual in mappy.fastx_read(infile, read_comment=False):
        reads[name] = seq
    return reads

def read_UMI(umi_file):
    UMIdict={}
    for line in open(umi_file):
        a=line.strip().split('\t')
        name_root=a[0].split('_')[0]
        if len(a)==2:
            UMI=a[1]
        else:
            UMI='_'
        if UMI not in UMIdict:
            UMIdict[UMI]=[]
        UMIdict[UMI].append(name_root)
    return UMIdict

def determine_consensus(UMI,reads,fastq_reads,temp_folder,subsample,medaka):
    '''Aligns and returns the consensus'''
    racon='racon'
    type1='best'
    corrected_consensus = ''
    fastq_count=0
    og_reads=reads
    temp_files=[]
    if len(reads)==1:
        for name,seq,qual,UMI in reads:
            main_name=name
            corrected_consensus=seq
    else:
        accuracies=[]
        out_Fq=temp_folder+'/'+UMI+'.fastq'
        temp_files.append(out_Fq)
        out=open(out_Fq,'w')

        indexes = np.random.choice(range(0, len(fastq_reads), 1),
                                   size=min(len(fastq_reads),subsample),
                                   replace=False)
        indexes=set(indexes)
        index=0
        for read in fastq_reads:
            if index in indexes:
                out.write('@'+read[0]+'\n'+read[1]+'\n+\n'+read[2]+'\n')
            index+=1
        out.close()

        poa_cons = temp_folder + '/'+UMI+'_poa_consensus.fasta'
        output_cons = temp_folder + '/'+UMI+'_corrected_consensus.fasta'
        overlap = temp_folder +'/'+UMI+'_overlaps.paf'
        temp_files.append(poa_cons)
        temp_files.append(output_cons)
        temp_files.append(overlap)
        overlap_fh=open(overlap,'w')

        names,max_coverage,repeats,qual,raw,before,after=[[],0,0,[],[],[],[]]

        for read,sequence,q,UMI in reads:
            first=sequence
            info=read.split('_')
            names.append(info[0])
            coverage=int(info[3])
            qual.append(float(info[1]))
            raw.append(int(info[2]))
            repeats+=int(info[3])
            before.append(int(info[4]))
            after.append(int(info[5].split('|')[0]))
            if coverage>=max_coverage:
                best=sequence
                best_name=read
                max_coverage=coverage

        ### running poa
        if len(reads)<=2:
            consensus_sequence=best
        else:
            sequences=[]
            mm_align = mappy.Aligner(seq=first, preset='map-ont')
            for name,sequence,q,UMI in reads:
                for hit in mm_align.map(sequence):
                    if hit.strand==1:
                        sequences.append(sequence)
                    elif hit.strand==-1:
                        sequences.append(mappy.revcomp(sequence))
                    if hit.mlen>0:
                        accuracies.append(1-(hit.NM/hit.mlen))

            if not sequences:
                print('no sequences','reads')
                consensus_sequence = best
            else:
                res = poa_aligner.msa(sequences, out_cons=True, out_msa=False)
                if not res.cons_seq:
                    print('abpoa not successful for',UMI,' - falling back on best R2C2 read',' '*20,end='\r')
                    consensus_sequence = best
                else:
                    consensus_sequence = res.cons_seq[0]
                    type1='poa'
        out_cons_file = open(poa_cons, 'w')
        out_cons_file.write('>Consensus\n' + consensus_sequence + '\n')
        out_cons_file.close()


        ### running racon

        final=poa_cons
        temp_cons_reads=read_fasta(poa_cons)
        mm_align = mappy.Aligner(seq=temp_cons_reads['Consensus'], preset='map-ont')
        for name,sequence,q in fastq_reads:
            for hit in mm_align.map(sequence):
                overlap_fh.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    name, str(len(sequence)), hit.q_st, hit.q_en,
                    hit.strand, 'Consensus', hit.ctg_len, hit.r_st,
                    hit.r_en, hit.mlen, hit.blen, hit.mapq))

        overlap_fh.close()
        os.system('%s -q 5 -t 1 --no-trimming\
                   %s %s %s >%s 2>./racon_messages.txt' \
                   %(racon,out_Fq, overlap, poa_cons, output_cons))

        final=output_cons
        reads = read_fasta(final)
        if len(reads)==0:
            fastq_count+=1
            print('racon not successful for',UMI,' - falling back on abpoa consensus',' '*20,end='\r')
            reads = read_fasta(poa_cons)

        else:
            type1='racon'
            UMIpath=temp_folder+'/'+UMI
            sub2draftSAM=UMIpath+'.sub2draft.sam'
            sub2draftBAM=UMIpath+'.sub2draft.bam'
            sub2draftBAMsorted=UMIpath+'.sub2draft.sorted.bam'
            sub2draftBAMsortedIndex=UMIpath+'.sub2draft.sorted.bam.bai'
            medakaHDF=UMIpath+'.medaka.hdf'
            medakaFasta=UMIpath+'.medaka.fasta'
            temp_files.append(sub2draftSAM)
            temp_files.append(sub2draftBAM)
            temp_files.append(sub2draftBAMsorted)
            temp_files.append(sub2draftBAMsortedIndex)
            temp_files.append(medakaHDF)
            temp_files.append(medakaFasta)

            if medaka:
                os.system('minimap2 -ax map-ont --secondary=no -t 1 %s %s >%s 2> ./minimap2.messages' % (final,out_Fq,sub2draftSAM))
                os.system('samtools view -b %s > %s'% (sub2draftSAM,sub2draftBAM))
                os.system('samtools sort %s > %s' % (sub2draftBAM,sub2draftBAMsorted))
                os.system('samtools index %s' % (sub2draftBAMsorted))
                os.system('medaka consensus %s %s --model r1041_e82_400bps_sup_v4.0.0 2>medaka_errors.txt > medaka_messages.txt' %(sub2draftBAMsorted,medakaHDF))
                os.system('medaka stitch %s %s %s 2>medaka_stitch.errors >medaka_stitch.txt' % (medakaHDF,final, medakaFasta ))
                final=medakaFasta
                reads=read_fasta(final)
                if len(reads)==0:
                    fastq_count+=1
                    print('medaka not successful for',UMI,' - falling back on racon consensus',' '*20,end='\r')
                    reads = read_fasta(output_cons)
                else:
                    type1='medaka'

        for read in reads:
            corrected_consensus = reads[read]
        main_name=best_name.split('_')[0]+'_'+str(round(np.average(qual),2))+'_'+str(int(np.average(raw)))+'_'+str(repeats)+'_'+str(int(np.average(before)))+'_'+str(int(np.average(after)))

    for temp_file in temp_files:
        if os.path.isfile(temp_file):
            os.system('rm '+temp_file)
    if accuracies:
        mean_acc=np.mean(accuracies)
    else:
        mean_acc=0
 #   print('finished UMI',UMI,mean_acc,' '*60,end='\r')
    return main_name,corrected_consensus,str(len(og_reads)),str(len(fastq_reads)),type1,mean_acc


def read_cons(fasta_file):
    ConsDict={}
    for name,seq,q in mappy.fastx_read(fasta_file):
        root=name.split('_')[0]
        ConsDict[root]=(name,seq,q)
    return ConsDict

def process_batch(subreads,reads,counter,processed,singles,subsample,medaka,threads,out,fasta_file,combined_UMI_length):
    results={}
    results1={}
    temp='./tmp'+fasta_file.split('.')[0]
    if not os.path.isdir(temp):
        os.system('mkdir '+temp)
    pool = mp.Pool(processes=threads)
    print(f'making consensus sequences in parallel for the next {len(reads)} UMIs', ' '*60)
    for temp_UMI,temp_reads in reads.items():
        temp_reads=list(temp_reads)
        UMIseq=temp_reads[0][3]
        if temp_UMI in subreads:
            temp_subreads=list(subreads[temp_UMI])
        else:
            temp_subreads=[]
        if len(temp_reads)>1 and len(UMIseq)==combined_UMI_length:
#            print('submitting UMI',temp_UMI,'with',len(temp_reads),'consensus reads and',len(temp_subreads),'total subreads',' '*20, end='\r')
            results[temp_UMI]=pool.apply_async(determine_consensus,[temp_UMI,temp_reads,temp_subreads,temp,subsample,medaka])
        else:
            counter=1
            for temp_read in temp_reads:
                results1[temp_UMI+'.'+str(counter)]=(temp_read[0],temp_read[1])
                counter+=1
    nameDict={}
    pool.close()
    pool.join()
    gc.collect()
    for temp_UMI,result in results.items():
        consName,consSeq,readNumber,subreadNumber,type1,accuracies = result.get()
        out.write('>%s_%s_%s_%s_%s_%s\n%s\n' % (consName,readNumber,subreadNumber,temp_UMI,str(round(accuracies,3)),type1,consSeq))
    for temp_UMI,result in results1.items():
        consName,consSeq = result
        out.write('>%s\n%s\n' % (consName,consSeq))
    return(processed+len(results),singles+len(results1))

def delegating_consensus_generation(fasta_file,umi_file,subread_file,subsample,medaka,combined_UMI_length,threads,batch,out):
    target=batch
    subsample=subsample
    pool = mp.Pool(processes=threads)
    print('reading UMI file')
    UMIdict=read_UMI(umi_file)
    totalUMIs=len(UMIdict)
    print(f'there are {len(UMIdict} UMIs in total')
    print('reading consensus reads')
    ConsDict=read_cons(fasta_file)
    counter,processed,singles,reads,subreads=0,0,0,{},{}
    print('partioning reads for multiprocessing')
    for line in open(subread_file):
        a=line.strip().split('\t')
        counter+=1
        name,seq,qual = a[2],a[3],a[4]
        UMInumber=a[0]
        UMI=a[1]
        if UMInumber not in subreads:
            if counter>target:
                print(f'\ncollected {counter} subreads covering {len(subreads)} UMIs', ' '*60)
                processed,singles=process_batch(subreads,reads,counter,processed,singles,subsample,medaka,threads,out,fasta_file,combined_UMI_length)
                print(f'{processed} reads merged from multi-read UMIs, {singles} reads not merged because incomplete or single-read UMI', ' '*60)
                subreads={}
                reads={}
                counter=0
            subreads[UMInumber]=set()
            reads[UMInumber]=set()
        rootNames=UMIdict[UMI]
        for rootName in rootNames:
            n,s,q = ConsDict[rootName]
            reads[UMInumber].add((n,s,q,UMI))
        subreads[UMInumber].add((name,seq,qual))

    processed,singles=process_batch(subreads,reads,counter,processed,singles,subsample,medaka,threads,out,fasta_file,combined_UMI_length)
    print(f'{processed} reads merged from multi-read UMIs, {singles} reads not merged because incomplete or single-read UMI', ' '*60)



