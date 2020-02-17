from sys import argv
import os
import re
import cut_adapt
import align
import paired_alignment
import generate_binned_midpoints
import count_interactions_per_binned_fragPairs
# integrated pipeline
## input two fastq files, output directory, reference genome, chromosome sizes, temp folder, adapter file
## output .hic files

## parameters and thresholds
BINSIZE = 5000
DISTANCETRHESHOLDS = 1000
MAPPABILITYTHRESHOLD = 1

USAGE = """
Usage: python process_hic.py <fastq1> <fastq2> <output directory> <bwa indexed reference genome> <chrome sizes> <temp directory> <adapters fasta file>

Process Hi-C raw data to 5kb resolusion contact file.
"""
if (len(sys.argv) != 8):
	  sys.stderr.write(USAGE)
	  sys.exit(1)

fastq1 = argv[1]
fastq2 = argv[2]
outDir = argv[3]
reference = argv[4]
chromsizes = argv[5]
temp = argv[6]
adapter = argv[7]

fileName1 = os.path.basename(re.search('(.*).fastq.*', fastq1).group(1))
fileName2 = os.path.basename(re.search('(.*).fastq.*', fastq2).group(1))

# step 1: adapter triming
print('step 1: adapter triming')
cut_adapt.main(fastq1, outDir, adapter)
cut_adapt.main(fastq2, outDir, adapter)
fastq1 = os.path.join(outDir, fileName1 + '.trimmed.fastq.gz')
fastq2 = os.path.join(outDir, fileName2 + '.trimmed.fastq.gz')

# step 2: alignment
print('step 2: alignment')
align.main(fastq1, outDir, reference, temp)
align.main(fastq2, outDir, reference, temp)
mapped1 = os.path.join(outDir, fileName1 + '.mapped')
mapped2 = os.path.join(outDir, fileName2 + '.mapped')

# step 3: paired alignment
print('step 3: paired alignment')
paired_alignment.main(mapped1, mapped2, outDir)
paired = os.path.join(outDir, 'paired_alignment_dedup')

# step 4: count binned interactions
print('step 4: count binned interactions')
binfile = paired + '.' + str(BINSIZE) + '.bins_midpoint'
intfile = os.path.join(outDir, 'interactions.txt')
mappability = os.path.join(outDir, 'mappability.txt')
generate_binned_midpoints.generateBinnedMids(BINSIZE, chromsizes, binfile)
count_interactions_per_binned_fragPairs.countInteractionAndEliminateNonMappables(paired, binfile, BINSIZE, MAPPABILITYTHRESHOLD, DISTANCETRHESHOLDS, intfile, mappability)
print('Finish!')
