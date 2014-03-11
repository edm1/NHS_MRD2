NHS_MRD2
========

## Instructions

Install the requirements with

`sudo python install.py`

Run the pipeline with

`python detect_targets_pipeline.py [raw_reads.fastq.gz] [output directory]`

Results files:

* target_results.txt - tab delimited results file
* NDN_clusters.fasta - clone N-D-N cluster sequences
* unmapped_reads.fasta - unmapped reads
* others

### Info

* Proportions in *target_results.txt* file are calculated using the total number of mapped reads.

### TODO

* Calculate and report consensus sequence for each cluster rather than using centroid sequence.
* Report the most frequently used V/J alignment for each cluster rather than the centroid alignment.
* At the moment CD-HIT is discarding any N-D-N sequences shorter than 10 bp. Need to change this.
* Make a log file.

* Implement follow-up stage

-------------------------------------

### Other

To process fastq using bowtie2:

1) Need to dereplicate the raw fastq.gz file

`python derep_fastq.py [in fastq.gz] [out fastq name]`

2) Need to use bowtie2 to search reads against (once for each index)

`bowtie2-align -p 4 --very-sensitive-local --reorder -x bowtie_indexes/J_w_phix_indexes -S J_out.sam -U [dereped fastq] `
`bowtie2-align -p 4 --very-sensitive-local --reorder -x bowtie_indexes/V_indexes -S V_out.sam  -U [dereped fastq]`

3) Need to process the SAM output files

`python process_sam.py [J.sam] [V.sam] [output_directory]`
