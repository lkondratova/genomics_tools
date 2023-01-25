import argparse
import pandas as pd
from collections import defaultdict
from gtfparse import read_gtf

parser = argparse.ArgumentParser()
parser.add_argument('--realtive positions', help = "File with columns 'gene', 'start', 'end'", required = True)
parser.add_argument('--gtf', help = "Self-explanatory", required = True)


args = parser.parse_args()
relative_positions = pd.read_csv(args.conserved_motifs, sep = '\t', low_memory = False) #gm_tfbs
gtf_all = read_gtf(args.gtf)

relative_positions = pd.read_csv('g_max_conserved_sites.txt', sep = '\t', low_memory = False) #gm_tfbs
gtf_all = read_gtf('Phaseolus_vulgaris.PhaVulg1_0.52.gtf')

#select genes only

gtf = gtf_all[gtf_all['feature'] == 'gene']

class gene_info:
    def __init__(self, gene, strand, start, end, chrom):
        self.gene = gene
        self.strand = strand
        self.start = start
        self.end = end
        self.chrom = chrom    

def add_gene_info(gene, strand, start, end, chrom):
    globals()[gene] = gene_info(gene, strand, start, end, chrom)
        
def gene_to_genomic_position(gene, gene_info, relative_start, relative_end, final_df):
    if gene_info.strand == '-':  #relative position is negative if upstream TSS and positive if downstream
        genomic_start = int(gene_info.start) + int(relative_start)
        if relative_end:
            genomic_end = int(gene_info.start) + int(relative_end)
        else:
            genomic_end = relative_end
    else:
        if relative_end:
            genomic_end = int(gene_info.end) - int(relative_start)
            genomic_start = int(gene_info.end) - int(relative_end)
        else:
            genomic_start = int(gene_info.end) - int(relative_start)
            genomic_end = relative_end        
    final_df. append([gene, gene_info.chrom, relative_start, relative_end, genomic_start, genomic_end])
    
    
    
gtf.apply(lambda x: add_gene_info(x['gene_id'], x['strand'], x['start'], x['end'], x['seqname']), axis = 1)  
final_df = []
relative_positions.apply(lambda x: gene_to_genomic_position(x['gene'], globals()[x['gene']], x['start'], x['end'], final_df), axis = 1)


df = pd.DataFrame(final_df)    
df.columns = ['gene', 'chromosome', 'relative_start', 'relative_end', 'genomic_start', 'genomic_end']   

df.to_csv("genomic_positions_out.txt", sep = '\t', header=True, index=False) 