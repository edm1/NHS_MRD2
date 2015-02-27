NHS MRD NGS Pipeline
====================
Version 1.5.3 - Pipeline for the analysis of minimal residual disease using NGS technologies.

(last updated 11/11/2013)

Installation
------------
The program *git* is required to track and apply updates that occur during development of the pipeline. It can be installed using `sudo apt-get install git`.

To retrieve the pipeline from github use `git clone https://github.com/edm1/NHS_MRD.git` (recommended) or download and unzip the [archive](https://github.com/edm1/NHS_MRD/archive/master.zip). The pipeline is designed to run on Ubuntu 12.04+ using python 2.7.4. It should work on any linux distribution but may require manual set up (see below).

#### Automatic set up (debain/ubuntu 12.04+)
The pipeline can be set up using the *setup_pipeline.py* script. For the intial set up run `sudo python setup.py install`.

#### Manual set up
To run the pipeline on a non-debain based linux distribution, the following python modules need to be installed:
- [matplotlib](http://matplotlib.org/)

You must also have docopt installed using python-pip:

```
sudo pip install docopt==0.6.1
```

And the following system packages must also be available in the system path:
- blastn (from the [BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) package)

Additionally, 'results', 'tmp', 'comparisons' and 'database' folders must be created and the [IgBLAST](http://www.ncbi.nlm.nih.gov/igblast/) databases must be downloaded. This can be done with the following commands:

```
rm -rf database tmp results comparisons
mkdir --parents database tmp results comparisons
wget -P ./database ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/database/human*
```

#### Updating the pipeline
If the pipeline was originally downloaded using the `git clone` command, it is possible to automatically update to the latest version using the command `python setup.py update` (this will overwite any changes you have made to the scripts themselves but will not affect any results files). Otherwise it will be necessary to re-download the latest files from github using the browser.

#### Updating the BLAST database
To re-download the latest BLAST database use the command `python setup.py updatedb`.

Usage
-----
### Pypy
The pipeline has been re-written to work using the *much* faster pypy interpretor (approximately 5 times faster). It is highly recommended you run the pipeline using `pypy` instead of `python`.

### 1. Detecting targets
The first stage in the analysis pipeline is to determine which clones are present in what quantities for each sample. This is executed using this command:

`pypy nhs_mrd.py detect <reads.fq.gz> <output_id>`

where \<reads.fq.gz> is location of the fastq file containing the reads in the direction from J -> V. The \<output_id> should be a unique identifer for this sample, it is used as a prefix for the results files.

#### Parameters
Additional parameters can be appended to the above command. To view all paramters and their default values use `python nhs_mrd.py --help`.

###### General

* `--v` or `--version` - Will show the current pipeline version
* `--top <INT>` - Number of top clusters to report in .tab and .detail results files (100)
* `--quiet` - Do not print progress to terminal.

###### Read filtering
* `--r <REGEX>` - Regular expression (regex) representing the motif used to identify J-read. This should be the reverse complement of the actual sequence. ("CCCCAG")
* `--m <INT>` - The maximum distance that the motif (--r) can be from the beginning of the read and still be considered a true motif (50)
* `--germ <REGEX>` - Regex that will match any germline sequences present in the sample (see appendix below)
* `--ab <INT>` - Number of bases to inculde after the successful motif (--r) match (125)
* `--bb <INT>` - Number of bases to inculde before the successful motif (--r) match. Note: this number should be smaller than the expected distance from the start of the read to --r motif (25)

###### Clustering

* `--min_dup <INT>` - After sequences are de-replicated, any groups containing less than this number of reads will be discarded.
* `--identity <INT>` - Should be an interger between 1 and 100. It represents the % min % indentity between two reads for it to be considered part of the same cluster.

###### Annotation
* `--t <INT>` - Number of BLAST hits to keep for V, D and J genes (3)
* `--e <DEC>` - E-threshold below which BLAST hits will be considered significant (0.05)

#### Input files
Input file for target detection is the raw fastq sequencing files. This can either be uncompressed or gzipped. Using default settings you should make sure there are at least 157 bases (from J end) with high PHRED scores.
* `--reads` are the J-end sequencing reads in J -> V direction

```
   |------------------V---------------|    |---D---|  |----J----|  
5'-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX-3'
                             <---------------------------reads--|

Ref:figure 1
```

#### Output files
Results from this stage of the pipeline are output into the *results* and *tmp* folders.

###### results folder
The *results* folder contains the main outputs from target detection.
* *prefix*.log - a log file showing: 
    - the parameters used
    - the number each read type present in sample (e.g. germline)
    - the number of clusters produced
* *prefix*.tab - summary of the top <N> (--top) clusters, including:
    - ID
    - number of reads
    - proportion of reads
    - top V/D/J genes
    - consensus sequence.
* *prefix*.detail - detailed results for each cluster, including:
    - differential V-usage
    - additional significant D/J hits
    - details of D/J/N rearrangments
    - shorthand rearrangement notation
* *prefix*.dat - data file required by the second stage of the pipeline.
* *prefix*.cons - fasta file containing the consensus sequences for **all** clusters. Sequence header shows:
    - the ID of the cluster's centroid
    - size of the cluster
* *prefix*_bar.pdf - 3D bar plot showing V vs. J gene usage

###### tmp folder
The *tmp* folder contains some additional files which may be of interest for troubleshooting:
* *prefix*_germline.fa - fasta file containing any detected germline sequences that were removed from the analysis
* *prefix*_Jreads.fa - fasta file containing trimmed V/J-end reads used in the analysis
* *prefix*.top_V/D/J.blast - xml file containing BLAST results from V/D/J segments
* *prefix*_motif_too_close.fa - fasta file containing reads where the conserved motif was too close to J-end.
* *prefix*_motif_too_far.fa - fasta file containing reads where the conserved motif was too far from J-end.
* *prefix*_no_motif.fa - fasta file containing reads that had no conserved motif.

### 2. Comparing samples
Once stage 1 has been carried out on two samples (i.e. a diagnosis and follow-up sample), the samples can be compared to find matching clones in each. To compare samples use:

```
pypy nhs_mrd.py compare <sample1.dat> <sample2.dat> <output id>
```

where \<sample1.dat> is the data file from clone detection in sample 1 (e.g. ./results/sample1/sample1.dat) and \<sample2.dat> is the equivalent for the other sample. Output ID is a unique identifier that will be used as a prefix for the comparison results.

#### Parameters
View all paramenters using `python nhs_mrd.py --help`.

* `--queries <INT>` - Number of top clusters to compare between each sample (10)
* `--matches <INT>` - Number of matches to report for each of the top clusters (3)
* `--quiet` - Do not print progress to terminal

#### Input files
Input files are the .dat files ( *prefix*.dat ) produced by the *detect_targets_pipeline.py* script and stored in the *results* folder.

#### Output files
Results from this second stage are stored in the *comparisons* folder. There is a single file produced:

* *prefix*.shared - For the top ten clusters in the diagnosis file, it shows the clipped sequence used as a query and clusters that have consensus sequences that match the query in both the diagnosis and follow-up samples. It then repeats the process with the samples switched over.

```
[sample1 cluster1] query sequence
    Matches in sample 1
        cluster_number mis-matches normalised_size actual_size proportion ID
        ...
        ...
        Total size of matches with mismatch of 0
        
    Matches in sample 2
        cluster_number mis-matches normalised_size actual_size proportion ID
        ...
        ...
        Total size of matches with mismatch of 0
        
    Fold-change between sample1 and sample 2
    
[sample1 cluster 2] ...
    ...
    ...
    
Ref:table 1
```

The columns in table 1 refer to:
* *cluster_number* - the cluster number for that cluster in the *.tab* results file. Cluster 1 is the largest cluster, cluster 2 is the 2nd largest and so on...
* *mis-matches* - the number of mismatches (actually edit distance) between the query sequence and the clusters consensus sequence
* *normalised_size* - that cluster's normalised sized as calculated using equation 1 (below)
* *actual_size* - the actual number of reads that made up that cluster
* *proportion* - the proportion of reads in the sample that made up that cluster
* *ID* - unique identifier for the cluster's centroid sequence. Same as the IDs in the original sequencing reads file

How the pipeline works
----------------------
### Detecting targets
The first script is split into 3 main stages.

#### (1/3) Read filtering
The sequences within the J-end fastq file ( *--reads* ) are searched for immunoglobulin motifs ( *--r* and *--germ* ). If a *--germ* motif is found the read is stored in the *prefix_germline.fa* file. If a *--r* motif is found and is within the correct distance from the start of the reads ( more than *--bb* and less than *--m* ) the J-end read is stored in *prefix_J.fa*. The sequence 125bp ( *--ab* ) upstream and 25bp ( *--bb* ) downstream of the motif ( *--r* ) is clipped out and passed to clustering.

```
   ----V----|        |----D----|   |----------------J-----------|  
5'-XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX-3'
      <-----------------------------ab (125)--|motif|--bb (25)-->
                 
Ref:figure 2
```

Reads than contain neither motif are logged but discarded from further analysis.

#### (2/3) Clustering
Clustering is carried out on spliced out sequences. Firstly, the reads are dereplicated so that any identical sequences are represented only once. For each dereplication group, the number of sequences that it contains and a common ID name is stored. Only groups with size greater than 1 are included in the analysis. The groups are sorted in order of decreasing size.

Clustering based on sequence similarity is also used to remove PCR and sequencing errors. The largest group is taken as the centroid sequence for a new cluster. This centroid is then compared to every other sequence that has not yet been clustered. If the distance between them is less that the amount allowed by *--identity* then the groups are merged into a cluster. Once all sequences have been compared to the centroid, the next largest group is selected as a new centroid.

Once all reads have been clustered, a consensus sequence is produced for each.

[1] *This essentially says that if there are only a small number of sequences in one group compared to the centroid then the distance between them is more likely to be caused by sequencing error and so should be merged.*

#### (3/3) Annotation
Annotation of reads is carried out using *blastn* from [NCBI BLAST+](http://www.ncbi.nlm.nih.gov/books/NBK1763/). BLAST parameters optimised for detecting V, D and J genes were taken from [Ye *et. al.* (2013)](http://www.ncbi.nlm.nih.gov/pubmed/23671333).

V, D and J genes are annotated from the consensus sequences generated during clustering. D gene hits are only assumed correct if they are not overlapping with V and J hits.

### Comparing samples
For each of the top clusters in the diagnosis sample a query sequence is produced. The query sequence is the consensus sequence that has been clipped before the V region. This is so that only the N1-D-N2-J region is used as a query. The query sequence is compared to each of the clusters in the diagnosis and follow-up samples. The diagnosis samples are searched as well as the follow-up sample in order to identify clones that been put in separate clusters due to V replacement causing the 5' end of the consensus sequence to vary.

Top matches (calculated as edit distances) to the query in the diagnosis and follow-up samples are reported. If more that one cluster has an exact match, the normalised sizes are summed for that sample. If exact matches are found in both samples then a fold-change difference is calculated for that query. 

The whole process is then repeated with diagnosis and follow-up samples swapped around.

**Note on normalisation:** normalisation of cluster sizes is done by transforming the actual cluster size by a function of the total number of reads sequenced. The function used is:

```
normalised size = actual size / (total number of reads sequenced x 10^-6)

Ref:equation 1
```

Considerations for the future
-----------------------------
#### Detecting other immunoglobulins and T-cells
When adapting the pipeline to detect antibodies other than IgG there are a few things to consider:
* The common motif used to detect true antibody reads during the initial filtering stage.
    - All IgGs have a single motif in common at their 3' end (CTGGGG). It will be necessary to find a motif that is common to the antibody of interest. If there is not a single motif then the pipeline will have to be adapted to work with a set of motifs.
    - A motif must be long enough to not be very likely to occur randomly in the sequence.
    - The position of the motif is also very important. Currently there is only a single primer being used for the 3' (J) end of the sequence. This means that the common motif should always be found in the same position. The position of the motif is used to confirm whether it is a true motif or occured randomly. The motif must be closer than *--m* bp and more than *-bb* from the end of the read. If multiple different primers are used for the J end or if multiple motifs are used then the motif distance will not be fixed and so more though will have to be put into the read filtering.
* The regular expression that detects germline sequences will have to be altered as these must be removed before analysis.
* Clustering currently relies on all of the post-filerting J end reads being both aligned and of the same lenght. They currently are because of the use of a common CTGGGG motif.

License
-------

Appendix
--------

Regular expression matching IgH germline sequences:

```
"CACGGTG[ATGCN]{22,24}AGAAACCCA|CACAGTC[ATGCN]{22,24}ACAAAAACA|CACACAG[ATGCN]{24,26}ACAAAAACC|CACACAG[ATGCN]{24,26}ACACAAACC|CACATTG[ATGCN]{22,24}ACAAAAACC|CACATTG[ATGCN]{20,22}ACAAAGAAC|CACAATG[ATGCN]{21,23}ACAAAAACC"
```

TODO
----
### High priority
* ~~Include portion of V in clustering. Annotate V using consensus sequence. Clip V before sample comparisons.~~
* Produce plot for sample comparisons (Clonal frequency in diagnosis vs followup)
* ~~Tidy up folder structure and have separate scripts for those that the user interfaces with.~~

### Medium priority
* ~~Find a way to save all results using pickle.dump() so that input to sample comparisons stage is simpler~~
* Use FastQC to produce a quality report for each file

### Low priority
* pickle.dump() the cluster and dereplication reference dictionaries. Allow extraction of all sequences for a single cluster, production of a HMM for the cluster.
* look into new method for heirarchical clustering, rather than the current centroid based clustering

Known bugs
----------
* Shorthand notation for V gene is currently incorrect. Not sure if it is a xml parsing error or a respresentation error.