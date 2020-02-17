from sys import argv
import os
import re

# step2
# alignment
## input trimmed_fastq, out_directory, indexed_reference_genome, temp_directory
## return sam_file, mapped_file
def main(fastq, outDir, reference, temp):
    fileName = os.path.basename(re.search('(.*).fastq.*', fastq).group(1))

    sai = os.path.join(outDir, fileName + '.sai')
    print('start bwa alignment')
    cmd = f'bwa aln -t 8 {reference} {fastq} > {sai}'
    print(cmd)
    os.system(cmd)
    print(f'sai file created for {fastq}')

    sam = os.path.join(outDir, fileName + '.sam')
    print('creating sam file')
    cmd = f'bwa samse {reference} {sai} {fastq} > {sam}'
    print(cmd)
    os.system(cmd)
    print(f'sam file created for {fastq}')

    mapped = os.path.join(outDir, fileName + '.mapped')
    awk = '($14=="X0:i:1" || $15=="X0:i:1") &&  $5>=30 && $12=="XT:A:U" {split($1,id,"#"); split($13,editDistField,":"); if (editDistField[1] == "NM" && editDistField[3] <= 3) {print id[1]\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$10}}'
    print('generating mapped reads')
    #OUTPUT MAPPED READS FROM SAM FILE
    cmd = f'cat {sam} | samtools view -S -F 4 - | awk \'{awk}\' | sort -T {temp} -k 1,1 > {mapped}'
    print(cmd)
    os.system(cmd)
    print(f'mapped file created for {fastq}')

if __name__ == '__main__':
    fastq = argv[1]
    outDir = argv[2]
    reference = argv[3]
    temp = argv[4]
    main(fastq, outDir, reference, temp)
