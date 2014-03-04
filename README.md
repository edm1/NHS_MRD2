NHS_MRD2
========

To process fastq using bowtie2:

1) Need to dereplicate the raw fastq.gz file

`python derep_fastq.py [in fastq.gz] [out fastq name]`

2) Need to use bowtie2 to search reads against (once for each index)

`bowtie2-align -p 4 --very-sensitive-local --reorder -x bowtie_indexes/J_w_phix_indexes -S J_out.sam -U [dereped fastq] `
`bowtie2-align -p 4 --very-sensitive-local --reorder -x bowtie_indexes/V_indexes -S V_out.sam  -U [dereped fastq]`

3) Need to process the SAM output files

`python process_sam.py [J.sam] [V.sam] [output_directory]`
