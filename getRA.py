import os
import argparse
from collections import defaultdict
from sys import argv
import operator
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from collections import OrderedDict
##The location and type of variation
def mutSeq(refSeq, readSeq):
    if len(refSeq) != len(readSeq):
        return 'not_equal *'
    mutPosNum = 0
    varList = []
    for i, (a, b) in enumerate(zip(refSeq, readSeq)):
        if a != b:
            mutPosNum += 1
            varList.append(f"{a}{i+1}{b}") 

    if mutPosNum == 0:
        return '0 WT'
    return f"{mutPosNum} {','.join(varList)}"

##Read sequences with fasta format
def readFastaFile(fastaFile):
    desc = ""
    seq = ""
    with open(fastaFile) as f:
        for line in f:
            if(line.startswith('>')):
                desc += line[1:].strip()
                continue
            seq += line.strip()
    return desc, seq
##Calculate the number of each variation based on ref seq
def read_variants(file_path, reference):
    variants = defaultdict(int)
    with open(file_path) as f:
        while True:
            line1 = f.readline().rstrip('\n')  #descriptor
            line2 = f.readline().rstrip('\n')  #Seq
            if not line1:
                break
            if len(line2) != len(reference):
                continue
            spt = line1.split('-')
            mutcount = int(spt[-1])
            mutinfo = mutSeq(reference, line2)
            num_mutations = int(mutinfo.split(' ')[0])
            if num_mutations > 10:
                continue
            variants[mutinfo] += mutcount
    return variants
def CountVariants(reference, before, after, output_path):
    before_variants = read_variants(before, reference)
    after_variants = read_variants(after, reference)
    with open(output_path, "w") as of:
        for mutinfo, before_count in before_variants.items():
            if before_count < 10:
                continue
            after_count = after_variants.get(mutinfo, 0)
            of.write(f"{mutinfo} {before_count} {after_count}\n")

    with open(output_path, "a") as of:
        for mutinfo, after_count in after_variants.items():
            if mutinfo not in before_variants:
                of.write(f"{mutinfo} 0 {after_count}\n")
## Calculate fitness based on mutant counts
def AnalysisRA(countFile, refSeq, raFile, cutoff):

    with open(countFile, 'r') as f:
        lines = f.readlines()

    wtActivity = 1.0
    sVarRAMap = defaultdict(float)
    pairRAMap = defaultdict(float)
    
    for line in lines:
        spt = line.split()
        if spt[0]=='0':
            wtActivity = float(spt[3])/float(spt[2])
            break
    with open(raFile, "w") as of:
        of.write(f'0\t{spt[2]}\t{spt[3]}\t{round(wtActivity, 3)}\n')
        # Single mutant
        m1List = [
            [
                '1',
                str(int(var[1:-1])),
                var[-1],
                str(before),
                str(after),
                str(round(after / before / wtActivity, 2))
            ]
            for line in lines if (spt := line.split())[0] == '1'
            if (before := float(spt[2])) >= cutoff
            and (after := float(spt[3]))
            and (var := spt[1])
            and (sVarRAMap[var] := after / before / wtActivity)
        ]
        for m1 in sorted(m1List, key=lambda x: int(x[1])):
            of.write('\t'.join(m1) + '\n')
        ## double mutants
        m2List = [
            [
                '2',
                str(int(var1[1:-1])),
                str(int(var2[1:-1])),
                var1[-1],
                var2[-1],
                str(before),
                str(after),
                str(round(sVarRAMap.get(var1, float('nan')), 2)),
                str(round(sVarRAMap.get(var2, float('nan')), 2)),
                str(round(after / before / wtActivity, 2))
            ]
            for line in lines if (spt := line.split())[0] == '2'
            if (before := float(spt[2])) >= cutoff
            and (after := float(spt[3]))
            and (varList := spt[1].split(','))
            and (var1, var2 := varList)
            and (pairRAMap[spt[1]] := after / before / wtActivity)
        ]
        for m2 in sorted(m2List, key=lambda x: (int(x[2]), int(x[1]))):
            of.write('\t'.join(m2) + '\n')

        m3List = [
            [
                spt[0],
                *[str(int(var[1:-1])) for var in varList],
                ''.join(var[-1] for var in varList),
                str(before),
                str(after),
                *[
                    str(round(sVarRAMap.get(var, float('nan')), 2))
                    for var in varList
                ],
                ','.join(
                    str(round(pairRAMap.get(f"{var1},{var2}", float('nan')), 2))
                    for i, var1 in enumerate(varList)
                    for var2 in varList[i + 1:]
                ),
                str(round(after / before / wtActivity, 2))
            ]
            for line in lines if (spt := line.split())[0] > '2'
            if (before := float(spt[2])) >= cutoff
            and (after := float(spt[3]))
            and (varList := spt[1].split(','))
        ]
        for m3 in sorted(m3List, key=lambda x: (int(x[3]), int(x[2]), int(x[1]), int(x[0]))):
            of.write('\t'.join(m3) + '\n')

def main(args):

    fasta_path = args.native_seq 
    file_before = args.before
    file_after = args.after
    outputfile = args.output_file
    min_reads = args.min_reads

    desc, ref = readFastaFile(fasta_path)
    output_dir_base = os.getcwd()
    output_path1 = os.path.join(output_dir_base, outputfile+".count")
    CountVariants(ref, file_before, file_after, output_path1)
    output_path2 = os.path.join(output_dir_base, outputfile+".ra")
    AnalysisRA(output_path1, ref, output_path2, min_reads)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument(
      "-n", "--native_seq", type=str,
      help="Path to native sequence(for RNA U should be T)"
    )
    parser.add_argument(
      "-b", "--before", type=str, default=None,
      help="""Path to the file containing the read counts before selection."""
    )
    parser.add_argument(
      "-a", "--after", type=str, default=None,
      help="""Path to the file containing the read counts after selection."""
    )
    parser.add_argument(
      "-min", "--min_reads", type=float, default=10.0,
      help="""Cutoff for the number of reads before selection."""
    )
    parser.add_argument(
      "-o", "--output_file", type=str, default=None,
      help="""Prefix of the file in which to store the output""",
    )
    args = parser.parse_args()

    main(args)


