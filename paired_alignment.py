from sys import argv
import os
import re

# step3
# paired alignment
## input R1_mapped, R2_mapped, out_directory
## return filtered, joined mapped file
def main(map1, map2, outDir):
    out = os.path.join(outDir, 'paired_alignment')
    out_dedup = out+'_dedup'

    # pairing alignment
    cmd = f'join {map1} {map2} > {out}'
    print(cmd)
    print('pairing alignment')
    os.system(cmd)

    # remove duplicates
    awk1 = '{s1 =\"+\"; s2 =\"+\"; if ($2 ==\"16\") {s1 =\"-\"}; if ($6 ==\"16\") {s2 =\"-\"}; print $3\"\\t\"$4\"\\t\"s1\"\\t\"$7\"\\t\"$8\"\\t\"s2\"\\t\"$5\"\\t\"$9\"\\t\"$1}'
    awk2 = '{print $9\"\\t\"$3\"\\t\"$1\"\\t\"$2\"\\t\"$7\"\\t\"$6\"\\t\"$4\"\\t\"$5\"\\t\"$8}'
    cmd = f'python ./sort_loci.py {out} 2 | awk \'{awk1}\' | sort -u -k 1,6 | awk \'{awk2}\' > {out_dedup}'
    print(cmd)
    print('removing duplicate')
    os.system(cmd)

if __name__ == '__main__':
    map1 = argv[1]
    map2 = argv[2]
    outDir = argv[3]
    main(map1, map2, outDir)
