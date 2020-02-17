# DNase I Hi-C Processing Pipeline
Pipeline for Hi-C raw data processing. Number of binned interactions will be counted in a txt file.

## Dependencies
1. python 3.7.3
2. cutadapt 2.6 with Python 3.7.3
3. bwa 0.7.17-r1188
4. samtools 1.9

## Installation
Include all the .py scripts in a same path.

## Usage
```python
python process_hic.py <fastq1> <fastq2> <output directory> <bwa indexed reference genome> <chrome sizes> <temp directory> <adapters file>
```
To change resolution (bin size), minimum distance and mappability threshold, please modify line 14-16 in `process_hic.py`.

* fastq1, fastq2:\
    Paired-end sequencing data. Can be zipped in .gz format.

* output directory:\
    A folder that contains all the outputs.

* bwa indexed reference genome:\
    Reference genome in fasta format. Need to be indexed under the same path:
    ```bash
    bwa index <ref.fa>
    ```

* chrome sizes:\
    A tab separated file recording sizes of chromatins in the reference genome, with its first column containing chromatin name and second column containing corresponding sizes. An example for chrome sizes file see `mm10.chrome.sizes_rmchr`.

* temp directory:\
    A folder storing temporary files.

* adapter file:\
    A file containing adapters in fasta format. An example for adapter file see `adapters.fa`.

## Output
* `<filename>.trimmed.fastq.gz`:\
  Raw data with adapters trimmed.
* `<filename>.trimmed.sai` & `<filename>.trimmed.sam`:\
  Alignment for each fastq file.
* `<filename>.trimmed.mapped`:\
  Mapped file for each fastq file.
* `paired_alignment` & `paired_alignment_dedup`:\
  Integrated paired alignment files for pair-end fastq data.
* `paired_alignment_dedup.<bin size>.bins_midpoint`:\
  Bin mid-points in reference genome.
* `interactions.txt`:\
  Tab-separated file recording interaction counts between pairs of bins, in the format of:\
  ```chr1       fragmid1  chr2      fragmid2  noofinteractions```
* `mappability.txt`:\
  File for 'fit-hic' analysis and has some unnecessary columns, in the format of:\
```chr       0       fragmid    number-of-interactions-this-fragment-participates-in       mappable-or-not(0/1)```