# BC1
Create consensus reads from UMI labeled long-reads

###Quickstart

```bash

sh setup.sh

python3 bc1.py -i 16S_C3POa.fasta -o 16S.BC1_consensus -s R2C2_Subreads.fastq -f -t 60 -u 5.38:50.ACAG.BDHVBDHVBDHV.AG,3.38:50.ACAG.BDHVBDHVBDHV.CG

```


```bash

  -i INPUT_FASTA, --input_fasta 
                        demultiplexed,trimmed, and reoriented R2C2 consensus reads as output by C3POa_postprocessing.py
  -o OUTPUT_FILE_ROOT, --output_file_root
                        root file name. Intermediate and final output files will be start with this
  -s SUBREAD_FILES, --subread_files
                        R2C2 subreads as output by C3POa. Can be multiple comma separated files
  -m, --medaka          
			if set, medaka will be run. required for highest accuracy but definitely slows things down. assumes data was produced with R10.4 pores at 400bp/s speed
  -f, --fuzzy           
			if set, reads with a single mismatch between their UMIs will be combined. By default, only reads with identical UMIs will be combined.
  -S STARSOLO_SAM_FILE, --STARsolo_sam_file
                        experimental setting.Probably best left alone for nowd bares can be used at any position
  -b SUBSAMPLE, --subsample
                        subreads will be subsampled to this number.
  -t THREADS, --threads
                        defines the number of threads the multiprocessing will use
  -u UMIPATTERNS, --UMIpatterns
                        UMI patterns separated by commas. An example would be '5.0:3.GACAG.NNNNNNNNNNNNNN.,3.0:3.CAC.NNNNNN.TTTT' that would indicate two UMIs. The first at the 5prime end of the read starting somewhere in the first 3 bases of the
                        read and is flanked on the left with 'GACAG' and is 14nt long. The second at the 3prime end of the read starting somewhere in the first 3 bases of the (reverse complemented) read and is flanked on the left with 'CAC' and on
                        the right with 'TTTT' and is 6nt long. IUPAC wild card bares can be used at any position

```



![Asset 8](https://github.com/christopher-vollmers/BC1/assets/28308271/3aaab974-07ec-4c08-868f-9c4887d60a0c)


