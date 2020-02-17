from sys import argv
import os
import re

# step1
# cutadapt
## input fastq, output_directory, adapter_fasta
## return trimmed fastq

#ADAPTOR SEQUENCE IS HARD CODED (BRIDGE ADAPTOR FOR DNASE HI-C)
#parameter -m to cutadapt is the minimum sequence to retain after adaptor is excised
#-m 20 means if the sequence is shorther than 20bp after trimming, we remove that read
#parameter -u trims a fixed number bases from 5 prime or 3 prime end of the read
#-u -25 corresponds to trimming 25 bases from 3 prime end


def main(inputFastq, outDir, adapter):
    fileName = os.path.basename(re.search('(.*).fastq.*', inputFastq).group(1))
    outputFastq = os.path.join(outDir, fileName+'.trimmed.fastq.gz')
    # outputReport = os.path.join(outDir, fileName +'.report')

    cmd = f'cutadapt -a file:{adapter} -m 20 -u -25 --cores=4 {inputFastq} -o {outputFastq}'
    print(cmd)
    os.system(cmd)

if __name__ == '__main__':
    #Input fastq file
    inputFastq = argv[1]
    #The folder to output results
    outDir = argv[2]
    adapter = argv[3]
    main(inputFastq, outDir, adapter)
