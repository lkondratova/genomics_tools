
"""
Compare two genome assemblies

Usage:
--ref1        genome 1 in fasta format (sorted)
--ref2        genome 2 in fasta format (sorted)
--prefix      prefix that should be added to chromosome names in genome 2 if chromosome names are different in two assemblies
              not required (for example chromosome 1 in genome 1 is 'Chr1' and in genome 2 '1', in this case --prefix Chr)
--addzero     adds zero in chromosome names 1-9 in genome 2 (for example, chromosome 1 in a genome 1 is '01' and in a genome 2 '1',
              --addzero parameter allows the program to recognize 1 and 01 as the same names), not required
--remzero    removes zero in chromosome names 1-9 in genome 2, not required
--ignore_names    do not apply prefix to records containing this name (for example --ignore_names scaffold will omit all scaffolds when applying prefix)
-- out        output file name without extension

Example 

python compare_two_assemblies.py --ref1 path/to/ref1
                                 --ref2 path/to/ref2
                                 --prefix Chr
                                 --addzero True
                                 --ignore_names scaffold
                                 
Output example:
ref1_chr ref1_chr_len ref2_chr  ref2_chr_len    identical
Chr01   51433939        1       52205531        False
Chr02   49670989        2       49040938        False
Chr03   53438756        3       52284309        False
Chr04   48048378        4       45960019        False
Chr05   40923498        5       40819286        False
Chr06   31236378        6       31977256        False
Chr07   40041001        7       51758522        False
Chr08   63048260        8       59662532        False
Chr09   38250102        9       37469608        False
"""

import argparse
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

parser = argparse.ArgumentParser()
parser.add_argument('--ref1', help = 'path to reference 1', required = True)
parser.add_argument('--ref2', help = 'path to reference 2', required = True)
parser.add_argument('--prefix', help = 'Prefix to change in reference 2 so that chromosome names match', required = False, default='')
parser.add_argument('--addzero', help = 'Add zero before chromosome numbers 1-9 in ref2 if absent in ref1', action = "store_true", default = False)
parser.add_argument('--remzero', help = 'Remove zero before chromosome numbers 1-9 in ref2 if absent in ref1', action = "store_true", default = False)
parser.add_argument('--ignore_names', help = 'Do not apply prefix, suffix, and zero rules to chromosemes containing this keyword (for example "scaffold")',
                   required = False, default = None)
parser.add_argument('--out', help='Output file name', required = False, default = 'output')

def change_names(ref1_names, ref2_names, prefix, addzero, remzero, ignore_names):
    chrom_dict = defaultdict(lambda:'') 
    ref2_names_t = [i for i in ref2_names]
    for i in range(0,9):
        if addzero:
            ref2_names_t[i] = '0' + ref2_names[i]
        if remzero:
            ref2_names_t[i] = ref2_names[1:] 
    for i in range(len(ref2_names_t)):
        name = ref2_names_t[i]
        if ignore_names:
            j = prefix + name if not ignore_names in name else name
        else:
            j = prefix + name
        chrom_dict[j] = ref2_names[i]
        #ref2_names_t.remove(i)
    ref1_unique = list(set(ref1_names) - set(chrom_dict.keys()))
    ref2_unique = list(set(chrom_dict.keys()) - set(ref1_names))
    return(chrom_dict, ref1_unique, ref2_unique)

def compare_chromosomes(ref1_names, ref1_dict, ref2_dict, chrom_dict, ref1_unique, ref2_unique, out):
    for i in ref1_names:
        if i in ref1_unique:
            out.append([i, len(ref1_dict[i]), 'None', 'None', 'False'])
        else:
            ref2_name = chrom_dict[i]
            out.append([i, len(ref1_dict[i]), ref2_name, len(ref2_dict[ref2_name]), 
                        'True' if ref1_dict[i]==ref2_dict[ref2_name]  else 'False'])
    for i in ref2_unique:
        out.append(['None', 'None', i, len(ref2_dict[i]), 'False'])
    return(out)

def default_d():
    return('')

args = parser.parse_args()
ref1_file = args.ref1
ref2_file = args.ref2
prefix = args.prefix
addzero = args.addzero
remzero = args.remzero
ignore_names = args.ignore_names

ref1 = SeqIO.parse(open(ref1_file),'fasta')
ref2 = SeqIO.parse(open(ref2_file),'fasta')
ref1_dict = defaultdict(default_d)
ref2_dict = defaultdict(default_d)

for fasta in ref1:
    ref1_dict[fasta.id.split(' ')[0]] = str(fasta.seq)
    
for fasta in ref2:
    ref2_dict[fasta.id.split(' ')[0]] = str(fasta.seq)


ref1_names = list(ref1_dict.keys())
ref2_names = list(ref2_dict.keys())
chrom_dict, ref1_unique, ref2_unique = change_names(ref1_names, ref2_names, prefix, addzero, remzero, ignore_names)
out = [['ref1_chr', 'ref1_chr_len', 'ref2_chr', 'ref2_chr_len', 'identical']]
out = compare_chromosomes(ref1_names, ref1_dict, ref2_dict, chrom_dict, ref1_unique, ref2_unique, out)
df = pd.DataFrame(out)
df.to_csv(args.out + '.tsv', sep = '\t', header=False, index=False)

