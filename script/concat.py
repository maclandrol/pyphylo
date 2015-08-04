#!/usr/bin/env python
import argparse
from glob import glob
import os
try:
    from Bio.SeqRecord import SeqRecord
    from Bio import AlignIO
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.Alphabet import generic_protein
except ImportError:
    raise ValueError('Requirement not met!, Biopython not found, exiting')


def concatenateSequences(geneorder, speclist, alignment, missing='-', sep='#', gp=0, alpha=generic_protein):
    """Concatenate alignment"""

    align_dict = {}
    for gene in geneorder:
        if gene in alignment.keys():
            cur_al = alignment[gene]
            al_len = cur_al.get_alignment_length()
            spec_dict = dict((x.id.split(sep)[1-gp], x) for x in cur_al)
            for spec in speclist :
                # default , put missing data
                adding = SeqRecord(Seq(missing*al_len, alpha))
                if(spec in spec_dict.keys()):
                    adding = spec_dict[spec]
                try :
                    align_dict[spec] += adding
                except (KeyError, ValueError):
                    align_dict[spec] = adding

    for spec in align_dict.keys():
        align_dict[spec].id = spec

    return MultipleSeqAlignment(align_dict.values())

if __name__ == '__main__':
    
    seqformat = ['fasta', 'maf', 'clustal', 'phylip-relaxed', 'nexus', 'stockholm']
    parser = argparse.ArgumentParser(description="concat - A script for sequence alignment concatenation.", epilog="See the whole program at https://github.com/maclandrol/pyphylo")

    parser.add_argument("-wd", '--aligndir', required=True, dest='workdir', help="Working directory containing every alignment")
    parser.add_argument("--species", "-s", type=argparse.FileType('r'), dest='species', help="Use only this species if provided")
    parser.add_argument("--order", type=argparse.FileType('r'), help="concatenation order - default : alphabetic order")

    parser.add_argument('-inf', "--inform", default='fasta', choices=seqformat, dest='inform', help="Alignment input format (each file should end by '.format'")
    parser.add_argument('-outf', "--outform", default='phylip-relaxed', choices=seqformat, dest='outform', help="Output file in which the concatenated alignment should be written")

    parser.add_argument('-o', "--output", default='concat', dest='output', help="Output file in which the concatenated alignment should be written")

    parser.add_argument("-gp", '--genepos', dest='gp', type=int, default=0, help="Gene position according to the separator in each header, 1 (post) or 0(pre)")
    parser.add_argument("--sep", default='#', help="Gene-specie separator : default '#'")


    args = parser.parse_args()


    align_glob = glob(args.workdir+"/*"+args.inform)
    aligndict = {}
    for path in align_glob:
        if(os.path.isfile(path)):
            c_basename = os.path.basename(path)
            if(c_basename):
                gene = c_basename.split('.')[0]
                aligndict[gene] = AlignIO.read(path, args.inform)

    # getting the geneorder 
    geneorder = []
    if(args.order):
        geneorder = [line.strip() for line in args.species if not line.startswith('#')]
    else:
        geneorder = sorted(aligndict.keys())


    # getting the specielist
    specielist = []
    if(args.species):
        specielist = [line.strip() for line in args.species if not line.startswith('#')]

    else:
        specielist = [x.id.split(args.sep)[1-args.gp] for align in aligndict.values() for x in align ]
        specielist = list(set(specielist))


    concat_align = concatenateSequences(geneorder, specielist, aligndict, sep=args.sep, gp=args.gp)
    AlignIO.write(concat_align, os.path.join(args.workdir, args.output), format=args.outform)