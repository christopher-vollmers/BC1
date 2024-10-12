[![Github release](https://img.shields.io/github/tag/christopher-vollmers/BC1.svg?label=Version)](https://github.com/christopher-vollmers/BC1/tags)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](http://perso.crans.org/besson/LICENSE.html)


# BC1
Create consensus reads from UMI labeled long-reads.
BC1 will take UMI-containing R2C2 consensus reads and their subreads as input to generate highly accurate R2C2+UMI reads. 
BC1 can also take regular ONT reads fastqs as input (-i). In that case, just omit the subreads flag (-s)

## Background ##

BC1 grew from a collection of tools for our work on single cell full-length cDNA sequencing using 10X in this [paper](https://doi.org/10.1186/s13059-022-02615-z)
Then we compiled it into a tool for our work on ultra-accurate amplicon sequencing in this [paper](https://doi.org/10.1093/pnasnexus/pgae336)

For both papers, we had sequencing reads that contained UMIs and we wanted to combine reads based on those UMIs into more accurate consensus reads. 
To be able to work with 10x fl-cDNA and Illumina style dual UMI amplicons libraries we designed a straightforward syntax for the parsing of all (most?) UMIs.

We designed BC1 to scale with ever increasing read numbers in long-read sequencing studies. 
While we use it mostly with ONT-based R2C2 data, it can be used on regular ONT reads and PacBio CCS reads. In both of those cases, subreads should not be submitted. In the ONT case because there are no subreads.
In the PacBio case because the subreads naming scheme hasn't been implemented in BC1 yet. Mostly though, PacBio subread files are enormous and therefore take forever to parse.

So, whether you have amplicon data, RNA-seq data, fl-cDNA data, or any other type of single read data (doesn't even have to be long-read data) that contains UMIs, BC1 provides you with a single command way of parsing and combining your reads. 

BC1 uses 

[abpoa](https://github.com/yangao07/abPOA) to make a preliminary multi sequence alignment and consensus read
[racon](https://github.com/isovic/racon) to polish that consensus read with subreads (if provided) or input reads

and optionally
[medaka](https://github.com/nanoporetech/medaka) to extra polish (slowly) the polished consensus read with subreads (if provided) or input reads
 
Running medaka to polish the output for highest accuracy only makes sense if the data is from ONT sequencers. For now, the pore/chemistry is hardcoded for R10.4/LSK114. It can be change in the code and we plan to ultimately have a flag to set it.
 

## Setup ##

You'll need minimap2 and samtools installed in your path. Everything else should be installed with the setup script. 

```bash

sh setup.sh

```

## Short example ##



Dual UMI amplicons that aren't trimmed and therefore have their UMI fat (38nt+) into the read. -f allows one mismatch in the UMIs :

```bash

python3 bc1.py -i 16S_C3POa.fasta -o 16S.BC1_consensus -p BC1_output/ -s R2C2_Subreads.fastq -f -t 60 -u 5.38:50.ACAG.BDHVBDHVBDHV.AG,3.38:50.ACAG.BDHVBDHVBDHV.CG

```

Single UMI fl-cDNA without subreads with a 5' UMI containing a sequence spacer and using the 3' sequence of the transcript as an extra "UMI" requiring perfect UMI matches

```bash

python3 bc1.py -i cDNA.trimmed.fasta -o cDNA.trimmed.bc1 -p BC1_output/ -t 30 -u 5.0:12.ACAG.NNNNNTGTTCTGATTNNNNN.TGGT,3.0:10.T.NNNNNNNNNNNNNNN.

```

Standard dual UMI Illumina IGH library allowing for one mismatch in the UMIs and running medaka for polishing

```bash

python3 bc1.py -i IGH.fasta -o IGH.bc1 -p BC1_output/ -u 5.0:12.GACAG.NNNNNNNNNNNNNN.,3.0:12.GACAG.NNNNNNNNNNNNNN. -s R2C2_subreads.fastq -f -t 60 -m

```

## Long Explanations ##

BC1 assumes you have a 

1) Consensus read file (in either fasta or fastq format)
2) Subread file (in fastq format)

Those can be produced by produced by either R2C2 on ONT sequencers or PacBio Iso-Seq or Kinnex

If you want to use medaka for polishing (-m), it also assumes that the data was produces with R10.4/LSK114 pore/chemistry on ONT sequencers. 

ran C3POa to process R2C2 data. It also assumes that That produces a consensus 

To use BC1, you need to know the sequences surrounding your UMIs. So take a look at the protocol. If you have the bases surrounding your UMI, this is how you structure the syntax. Individual information is separated by periods.

1) What end of the read to look - 5' or 3'with 5' being at the start. For 3', the read is reverse complemented before looking
2) What position to start looking:What position to stop looking for the UMI.
3) Bases flanking your UMI on the left
4) Bases of your UMI (Accepts N, A,T,C,G, and IUPAC base symbols)
5) Bases flanking your UMI on the right

So if your read looks like this

ATCGTCTTTNNNNANNNNCCCCC

your syntax could look like this

5.4:10.TTT.NNNNANNNN.CCC

if you have more than one UMI, just add another syntaxs separated by a commas


## Run options ##


```bash

  -i, --input_fasta 
                        demultiplexed,trimmed, and reoriented R2C2 consensus reads as output by C3POa_postprocessing.py
  -o, --output_file_root
                        root file name. Intermediate and final output files will be start with this
  -p, --output_path
                        output directory. Intermediate and final output files will be in this directory
  -s, --subread_files
                        R2C2 subreads as output by C3POa. Can be multiple comma separated files
  -m, --medaka          
			if set, medaka will be run. required for highest accuracy but definitely slows things down. assumes data was produced with R10.4 pores at 400bp/s speed
  -f, --fuzzy           
			if set, reads with a single mismatch between their UMIs will be combined. By default, only reads with identical UMIs will be combined.
  -b, --subsample
                        subreads will be subsampled to this number.
  -t, --threads
                        defines the number of threads the multiprocessing will use
  -u, --UMIpatterns
                        UMI patterns separated by commas. An example would be '5.0:3.GACAG.NNNNNNNNNNNNNN.,3.0:3.CAC.NNNNNN.TTTT' that would indicate two UMIs. The first at the 5prime end of the read starting somewhere in the first 3 bases of the
                        read and is flanked on the left with 'GACAG' and is 14nt long. The second at the 3prime end of the read starting somewhere in the first 3 bases of the (reverse complemented) read and is flanked on the left with 'CAC' and on
                        the right with 'TTTT' and is 6nt long. IUPAC wild card bares can be used at any position
  -r, --resume          
                        if set, bc1 will try to resume from previous output. Will only consider past results that were generate with the same exact command settings.


```



![Asset 8](https://github.com/christopher-vollmers/BC1/assets/28308271/3aaab974-07ec-4c08-868f-9c4887d60a0c)


