from Bio.Seq import Seq
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_protein

import sys, os

def gene_type(gene):
    return gene

def get_color(gene):
    if(gene.startswith('trn')):
        return 'yellow'
    elif(gene.startswith('orf')):
        return 'pink'
    elif(gene[0:3] in ['rpl', 'rrn', 'rnl', 'rns', 'rps']):
        return 'green'
    else:
        return 'blue'


def parse_master_file(masterfile):

    genes = {}
    genomic_title = ''
    genomic_seq = ''
    gene_list = []

    with open(masterfile) as IN :
        glist = False
        end = False
        for line in IN:

            if 'end' in line.strip():
                end = True

            elif glist and (not end) and line.strip() != ';;':
                gene_list.extend(line.lstrip(';;').strip().split())

            elif(line.startswith(';;') and "List of genes" in line):
                glist = True

            if (line.startswith('>')):
                genomic_title = line.strip()
                break

        curpos = 1
        strand = 1

        for line in IN:
            line = line.rstrip()
            if(line.startswith(';')):
                line = line.split('G-')[1]
                strand = 1 if '==>' in line else -1
                #print line.split()[:3]
                g, _ , state  = line.split()[:3]
                if g not in genes.keys():
                    genes[g] = SeqFeature(FeatureLocation(curpos, curpos), id=g, type=gene_type(g), strand=0)

                genes[g].strand = strand
                if(state == 'start'):
                    genes[g].location._start = curpos
                elif(state == 'end'):
                    genes[g].location._end = curpos

            else :
                curpos, seq = line.split()
                genomic_seq+=seq.strip()
                curpos = int(curpos) + len(seq.strip())


    genomerecord = SeqRecord(Seq(genomic_seq, generic_protein), id=genomic_title.split('_')[0], description=genomic_title)

    return genes, genomerecord
