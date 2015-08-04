#!/usr/bin/env python
from Bio.Seq import Seq
from Bio import SeqIO, AlignIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from reportlab.lib.pagesizes import A4, inch, landscape
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph
from reportlab.lib.styles import getSampleStyleSheet

from Bio.Graphics import GenomeDiagram
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_protein, IUPAC
from Bio import Entrez #Taxonomy lookup

from prettytable import PrettyTable, ALL
from Bio.Blast import NCBIWWW, NCBIXML
import numpy as np #Array and matrix sums
from scipy import percentile,mean,std
#import subprocess, threading #Background process class
from Bio.Align import MultipleSeqAlignment #Create an alignment (for MUSCLE)
from ete2 import PhyloNode
import time #For waiting between sequence downloads
import argparse #For command line arguments
import sys, os, re
import copy #Getting subtrees
import sqlite3
import random, shutil
import collections

from applications import *
from phydata.mtgenes import revmtgenes
from masterfileparser import parse_master_file, get_color

# Setting Entrez parameters
Entrez.email = email
Entrez.tool = "pyPhylo_LocalScript"
maxCheck = 2
ALIGN_METHODS = {'muscle':['muscle'], 'mafft': ['mafft'], 'hmmmuscle' : ['hmmmuscle'],\
                'clustalo': ['clustalo'], 'prank' : ['prank'],\
                 'fsa':['fsa'], 'fast': ['muscle', 'mafft', 'clustalo'],\
                 'multi' : ['muscle', 'mafft', 'clustalo', 'fsa'],\
                 'all': ['muscle', 'mafft', 'clustalo', 'prank', 'fsa']}


def regexp(expr, item):
    """Regexp for search in sqlite database"""
    reg = re.compile(expr)
    return reg.search(item) is not None

def taxonIDLookup(taxonID):
    """Lookup a taxon ID (an integer) in the NCBI taxonomy.
    Returns (Species_name, (taxonomy_genus, taxonomy_family, etc.))
    Will likely throw 'server errors' until intenal timeout is reached if given anything else."""
    finished = 0
    while finished <= maxCheck:
        try:
            handleDownload = Entrez.efetch(db="taxonomy", id=taxonID, retmode="xml")
            resultsDownload = Entrez.read(handleDownload)
            handleDownload.close()
            finished = maxCheck + 1
        except:
            if finished == 0:
                print("!!!Server error - retrying...")
                finished += 1
                time.sleep(2)
            elif finished == maxCheck:
                print("!!!!!!Unreachable. Returning nothing.")
                return(tuple())
            else:
                finished += 1
                time.sleep(2)
    scientificName = resultsDownload[0]['ScientificName']
    lineage = resultsDownload[0]['Lineage'].split("; ")
    lineage.reverse()
    lineage = tuple(lineage)
    taxId = resultsDownload[0]['TaxId']
    mitoCode = resultsDownload[0]['MitoGeneticCode']['MGCName']
    return(scientificName, lineage, taxId, mitoCode)


def commonLookup(spName):
    """Lookup a species according to its common (English?) name.
    Returns (Species_name, (taxonomy_genus, taxonomy_family, etc.))
    """
    finished = 0
    while finished <= maxCheck:
        try:
            handleSearch = Entrez.esearch(db="taxonomy", term=spName)
            resultsSearch = Entrez.read(handleSearch)
            handleSearch.close()
            finished = maxCheck + 1
        except:
            if finished == 0:
                print("!!!Server error checking " + spName + " - retrying...")
                finished += 1
                time.sleep(2)
            elif finished == maxCheck:
                print( "!!!!!!Unreachable. Returning nothing.")
                return(tuple())
            else:
                finished += 1
                time.sleep(2)
    if resultsSearch['IdList']:
        return taxonIDLookup(resultsSearch['IdList'][0])
    else:
        return(tuple())

def cladeSpecies(cladeName):
    """Finds information about a clade.
    Returns [(species_a, (genus_a, family_a, etc_a)), (species_b, (genus_b, family_b, etc_b), cladeID, mitcohondrial_code]"""
    searchTerm = cladeName + '[subtree] AND species[rank]'
    finished = 0
    while finished <= maxCheck:
        try:
            handleSearch = Entrez.esearch(db="taxonomy", term=searchTerm)
            resultsSearch = Entrez.read(handleSearch)
            finished = maxCheck + 1
        except:
            if finished == 0:
                print("!!!Server error checking" + cladeName+ " - retrying...")
                finished += 1
                time.sleep(2)
            elif finished == maxCheck:
                print ("!!!!!!Unreachable. Returning nothing.")
                return()
            else:
                finished += 1
                time.sleep(2)
    if resultsSearch['IdList']:
        output = []
        for spId in resultsSearch['IdList']:
            output.append(taxonIDLookup(spId))
        return output
    else:
        return()

def findLineage(spName):
    """Finds lineage (taxonomy) of a species.
    Returns [genus, family, ...]"""
    try:
        handleSpName = Entrez.esearch(db="taxonomy", term=spName)
        resultsSpName = Entrez.read(handleSpName)
        handleSpName.close()
        handleID = Entrez.efetch(db="taxonomy", id=resultsSpName['IdList'], retmode="xml")
        resultsSpID = Entrez.read(handleID)
        lineage = resultsSpID[0]["Lineage"].split("; ")
        lineage.append(spName)
        lineage.reverse()
        return lineage
    except:
        return ()

def isInt(value):
    """ Verify if a string
    can be casted to an int"""
    try:
        int(value)
        return True
    except ValueError:
        return False

def isFloat(value):
    """ Verify if a string
    can be casted to an int"""
    try:
        float(value)
        return True
    except ValueError:
        return False

def check_specie_list(specielist,  taxonID=False, get_clade=False):
    """ Verify and return the true name
    of each specie in your specielist """

    print( """
    --------------------------------------------
    Verify if your species are genuine species
    --------------------------------------------

    """)
    curated_spec_list = []
    for spec in specielist:
        spechecked = None
        if(taxonID and isInt(spec)):
            spechecked = taxonIDLookup(spec)
        else :
            spechecked = commonLookup(spec)
            if not spechecked and get_clade:
                time.sleep(1)
                spechecked = cladeSpecies(spec.split()[0])
        if type(spechecked) == list:
            for cspec in spechecked:
                curated_spec_list.append(cspec[0])
        elif spechecked:
            curated_spec_list.append(spechecked[0])
        else:
            print (spec + " was not added")
    return curated_spec_list


def blast_sequence(seqrecord, workdir=""):
    """ check some orf origin"""
    if(workdir):
        workdir+="/"
    handle = NCBIWWW.qblast("blastp", "nr", seqrecord.format("fasta"))
    filename = workdir+serqrecord.id+"xml"
    file = open(filename, "w")
    file.write(handle.read())
    file.close()
    hande.close()


def annote_genes(record, genes, genelist, format="circular", wd="Data/outpdf"):
    """ output an annotation pdf"""

    if not os.path.exists(wd):
        os.mkdir(wd)

    gd_diagram = GenomeDiagram.Diagram(record.id)
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()

    for gene in genelist :
        if len(gd_feature_set) % 2 == 0:
            color = "light"+get_color(gene)
        else :
            color = get_color(gene)
        gd_feature_set.add_feature(genes[gene], sigil="ARROW", arrowshaft_height=0.7,
                                   color=colors.__getattribute__(color), label=True,
                                   label_size = 6, label_angle=0)

    gd_diagram.draw(format="circular", circular=True, pagesize=(20*cm,20*cm),
                    start=0, end=len(record), circle_core= 0.5)
    gd_diagram.write(wd+'/annotation.pdf', "PDF")


class pData():
    """A simple class to save and load state and data"""
    def __init__(self, pyphy):
        self.genelist = pyphy.genelist
        self.speclist = pyphy.speclist
        self.gene2spec = pyphy.gene2spec
        self.spec2genes = pyphy.spec2genes
        self.editmethod = pyphy.editmethod
        self.alignmethod = pyphy.alignmethod
        self.alignment = pyphy.alignment
        self.concat = pyphy.concat
        self.concatpart = pyphy.concatpart
        self.edited = pyphy.edited
        self.oridata = (pyphy.oriseq, pyphy.orispec)
        self.tree = pyphy.tree
        self.order = pyphy.order
        self.state = pyphy.state
        self.concatfile = ""
        if hasattr(pyphy, 'concatfile'):
            self.concatfile = pyphy.concatfile

class pyPhylo(object):

    def __init__(self, wd, specielist, genbank, corrdb, oridata, complete=True, \
                        editmethod='both', alignmethod='muscle', orthomcl=True, \
                        prompt=True, concat_order=[]):
        """Initialisation of the pyphylo object"""
        self.spec2genes = collections.defaultdict(list)
        self.gene2spec = collections.defaultdict(dict)
        self.workdir = wd
        self.aligndir = os.path.join(wd, aligndir)
        self.fastadir = os.path.join(wd, fastadir)
        self.orthoinput = os.path.join(wd, ortho_indir)
        self.orthooutput = os.path.join(wd, ortho_outdir)
        self.outpdf = os.path.join(wd, outputpdf)
        self.concatfile = concatfile
        self.outtree = os.path.join(wd, outputtree)
        self.datadump = os.path.join(wd, datadumping)
        self.order = concat_order
        self.state = 0
        # create directories if not exist
        for d in [self.workdir, self.fastadir, self.orthoinput, self.orthooutput, self.outpdf, self.outtree, self.datadump]:
            if not os.path.exists(d):
                os.mkdir(d)

        if( not (prompt and self.ask_prompt(choices=['load'], init=True))):
            self.oriseq, self.orispec = oridata
            self.alignmethod = alignmethod
            self.editmethod = editmethod
            self.alignment = collections.defaultdict(dict)
            self.edited = collections.defaultdict(dict)
            self.concat = None
            self.concatpart = None
            self.tree = None
            all_seq = SeqIO.index_db(os.path.join(wd,dbdir,seqdb), genbank, format="genbank")
            self.speclist = specielist
            # add the original data to the dict list
            for g, rec in self.oriseq.items():
                if orthomcl or 'orf' not in g:
                    if g.lower() in revmtgenes.keys():
                        g = revmtgenes[g.lower()]
                    rec.name = g
                    self.spec2genes[self.orispec].append(rec)
                    rec.id = rec.name + "#"+"_".join(self.orispec.split())
                    print ("Added : %s" % rec.description)
                    rec.description = ""
                    self.gene2spec[g][self.orispec] = rec

            def motif_in_seq(motif, slist, product):
                return sum([1 for x in slist if (motif in x.lower()) ]) > 0 and \
                                    'hypothetical' in  " ".join(product).lower()

            def can_add_seq(seq, spec, gene, g2s):
                """ verify if we should add this new sequence """
                try:
                    curr_seq = g2s[gene][sgene]
                    if(len(seq)>len(curr_seq)):
                        return True
                    return False
                except :
                    return True

            # fetch gene and protein from the database and add them to the list
            with sqlite3.connect(corrdb) as conn:
                conn.create_function("REGEXP", 2, regexp)
                holder = "|".join(self.speclist)
                query = "SELECT seqid, organism FROM genbank WHERE organism REGEXP ?"
                if(complete):
                    query += " AND complete_genome = 1"
                cursor = conn.execute(query, [holder])
                counter = 0
                #print("Au debut "+ str(len(self.speclist)))
                for row in cursor:
                    s_id, spec = row[0], row[1]
                    spec = spec.replace('.', '').replace('-', '_')
                    curr_seq = all_seq[s_id]
                    for f in curr_seq.features:
                        # this should filter hypothetical protein that we do not want
                        if f.type.upper()=='CDS' and 'gene' in f.qualifiers.keys() and not motif_in_seq('orf', f.qualifiers['gene'], f.qualifiers['product']) :
                            sgene = f.qualifiers['gene'][0].lower()
                            if(sgene in revmtgenes.keys()):
                                sgene =  revmtgenes[sgene]

                            prot_id = sgene+"#"+"_".join(spec.split()) # f.qualifiers['protein_id'][0]
                            descr = "%s %s; %s"%(spec, "-".join(f.qualifiers['gene']), prot_id)
                            seq = None
                            if 'translation' in f.qualifiers.keys():
                                seq = Seq(f.qualifiers['translation'][0], IUPAC.protein)
                            else :
                                try:
                                    table = int(f.qualifiers['transl_table'][0])
                                    seq = f.extract(curr_seq.seq).translate(table=table, to_stop=True)
                                except :
                                    pass

                            seq = seq.strip('*')
                            if (seq and can_add_seq(seq, spec, sgene, self.gene2spec)):

                                print ("Added : %s" %descr)
                                rec = SeqRecord(seq, id=prot_id, name=sgene, description="")
                                intercept = False
                                for i, r in enumerate(self.spec2genes[spec]):
                                    if(r.id == rec.id):
                                        intercept = True
                                        self.spec2genes[spec][i] = rec
                                if not intercept:
                                    self.spec2genes[spec].append(rec)
                                self.gene2spec[sgene][spec] = rec


            error, binaries = check_binaries([orthomclpipe])
            if not error:
                self.orthomcl()
            self.speclist = self.spec2genes.keys()
            self.genelist = self.gene2spec.keys()
            self.state = 1
            #print ("A la fin " + str(len(specielist)))

        self.ask_prompt()

    def init_from_pk(self, pk):
        """initalisation from pickle data"""
        self.genelist = pk.genelist
        self.speclist = pk.speclist
        self.gene2spec = pk.gene2spec
        self.spec2genes = pk.spec2genes
        self.alignment = pk.alignment
        self.concat = pk.concat
        self.concatpart = pk.concatpart
        self.edited =  pk.edited
        self.oriseq, self.orispec = pk.oridata
        self.tree = pk.tree
        self.order = pk.order
        self.editmethod = pk.editmethod
        self.alignmethod = pk.alignmethod
        self.state = pk.state
        self.concatfile = ""
        if hasattr(pk, 'concatfile'):
            self.concatfile = pk.concatfile

    def orthomcl(self, confile=None):
        """
        run orthomcl
        """
        print("""
            --------------
            Run orthomcl
            --------------

            """)
        # remove everything from current indir
        self.purge_directory(self.orthoinput)

        clusters = run_orthomcl(self.spec2genes, indir=self.orthoinput,\
                outdir=self.orthooutput, confile=confile)

        for groups in clusters:
            gs = [x.split('#') for x in groups]
            genes,species = zip(*gs)
            if(len(set(genes)) > 1):
                predo_genes = collections.Counter()
                for g in genes:
                    predo_genes[g] = len(self.gene2spec[g])
                freq_g = predo_genes.most_common(1)[0][0]

                for i,g in enumerate(genes):
                    s = species[i].replace('_', ' ')
                    if (g != freq_g and s not in self.gene2spec[freq_g].keys()):
                        can_add = True
                        for sp in self.gene2spec[g].keys():
                            if(sp in self.gene2spec[freq_g].keys()):
                                can_add = False
                        if(can_add):
                            print(g + "  to  " + freq_g)
                            print( "species : " + str(self.gene2spec[g].keys()))
                            self.gene2spec[freq_g].update(self.gene2spec.pop(g))


    def ask_prompt(self, choices=None, init=False):
        """
        show prompt
        """
        if not choices:
            choices = ['edit', 'report', 'save', 'align', 'refine', 'reset', 'load', 'concat','distR', 'fusgene']
        actions = {'edit': self.edit_data, 'save': self.dump_data, \
                    'report' : self.report_data, 'load': self.load_data, \
                    'align': self.perform_align, 'refine' : self.refine_align, \
                    'concat': self.concat_align,  'distR' : self.distR, \
                    'fusgene': self.fusgene, 'reset': self.reset_data}

        print ("""
            --------------------------------------------------------------------
            Choose one of the following or press enter to skip:

            %s

            --------------------------------------------------------------------

            """%(' -- '.join(choices)))

        action = raw_input('waiting for your action : ')
        action = action.strip()
        while (action and action not in choices):
            action = raw_input('Incorrect choice, Enter action : ')
            action = action.strip()

        if(action):
            actions[action]()
            if(init):
                return True
            else:
                self.ask_prompt()

        else: return False


    def fusgene(self, gene1=None, gene2=None, order=None):
        """ Function to make a fusion of gene"""
        if gene1 is None:
            gene1 = raw_input("Enter first gene : ")
        if gene2 is None:
            gene2 = raw_input("Enter second gene : ")
        if (gene1 not in self.genelist) or (gene2 not in self.genelist):
            print("Genes not found, returning to main menu...\n")
        else:
            if self.state < 2 :
                self.alignment['muscle'] = self.alignSequences(method='muscle', \
                                                outdir=None, timeout=99999999, verbose=False)
            alignid = self.alignment.keys()[0]
            if not order:
                order = sorted([gene1, gene2])
            if(len(order) > 2):
                print("Order argument is bigger than two, try again\n")
                return
            else:
                genefusion, fusionpart = self.concatenateSequences(order, self.speclist, self.alignment[alignid], missing='-', alpha=generic_protein)
                if genefusion:
                    self.genelist = list(set(self.genelist) - {gene1, gene2})
                    new_gene = "_".join(order)
                    self.genelist.append(new_gene)
                    # print self.spec2genes.keys()
                    for seq in genefusion:
                        spec = seq.id.replace('_', ' ')
                        seq.name = new_gene
                        seq.id = '%s#%s'%(new_gene, spec)
                        uniq_char = "".join(set(str(seq.seq)))
                        if uniq_char != '-':
                            self.gene2spec[new_gene][spec] = seq
                            self.spec2genes[spec].append(seq)


    def distR(self, rate_thresh=-1, align=None):
        """load rate file and use it to filter gene"""
        print("This is experimental\n")

        if self.state < 2 :
            self.alignment['muscle'] = self.alignSequences(method='muscle', \
                                            outdir=None, timeout=99999999, verbose=False)
            align = 'muscle'

        if not align:
            align = self.selection_alignment(self.alignment, filename=filename)

        if not(align):
            print("Alignment not provided, nor selected, returning to menu...\n")
            return

        distdir = os.path.join(self.aligndir, align)
        executeCMD("script/distmat.sh %s"%(distdir))
        distFile = os.path.join(self.aligndir, align, "tree.rates")
        try:
            with open(distFile, 'r') as DISTFILE:
                for line in DISTFILE:
                    print(line.replace(".fasta", "\t"))
                print("\n")
            retry = True
            while(retry):
                rate_thresh = raw_input("Enter your required threshold : ")
                if not isFloat(rate_thresh):
                    retry = raw_input("Error : rate should be a float, retry (y/n) ? ")
                    retry = retry.lower().startswith('y')
                else:
                    retry = False
                    rate_thresh = float(rate_thresh)
            if(rate_thresh > 0):
                genes, rates, outliers = load_rate_file(distFile, rate_thresh)
                print("This genes will be discarded for being outliers : \n")
                print outliers
                outgenes = zip(genes[outliers], rates[outliers])
                for g, r in outgenes:
                    print("%s : %f\n"%(g, r))
                    confirm = raw_input("Please confirm y/n : ")
                    if(confirm.lower().startswith('y')):
                        self.genelist = list(set(self.genelist) - set(genes[outliers]))
        except :
            print("Error, (check your distFile), returning to menu\n")


    def reset_data(self):
        """reset data to global values"""
        self.speclist = self.spec2genes.keys()
        self.genelist = self.gene2spec.keys()

    def report_data(self):
        """report data into a matrice of specie X gene
        """
        self.report_sequence(show=True)


    def load_data(self):
        """
        load current state to continue previous analysis
        """
        try :
            import cPickle as pickle
        except ImporError:
            import pickle

        print( """
                ---------------------------
                Loading mode
                ---------------------------
            """)

        filename = None
        inp = raw_input("Enter filename (not complete path) or press enter : ")
        if inp:
            filename = inp
        else :
            filename ="pyphylo.pkl"

        with open(os.path.join(self.datadump, filename), 'rb') as output:
            pk = pickle.load(output)
            self.init_from_pk(pk)


    def dump_data(self):
        """
        dump current state to file and save
        """
        try :
            import cPickle as pickle
        except ImporError:
            import pickle

        print( """
                ---------------------------
                Saving mode
                ---------------------------
            """)

        # save state as a pickle to reload
        filename = None
        inp = raw_input("Enter filename (not complete path) or press enter : ")
        if inp:
            filename = inp
        else :
            filename ="pyphylo.pkl"


        pk = pData(self)

        with open(os.path.join(self.datadump, filename), 'wb') as output:
            pickle.dump(pk, output, pickle.HIGHEST_PROTOCOL)

        #save gene list and save alignment
        self.save_genes()
        self.save_alignment()
        self.save_concat()


    def save_genes(self, format="fasta"):
        """Save fastafile in a directory"""
        self.purge_directory(self.fastadir)
        if self.genelist:
            for gene in self.genelist:
                SeqIO.write(self.gene2spec[gene].values(), os.path.join(self.fastadir,gene+".fasta"), format=format)


    def save_alignment(self, format="fasta"):
        """Save alignment in a directory"""
        # save current genelist also
        if self.alignment:
            for method in self.alignment.keys():
                methoddir = os.path.join(self.aligndir,method)
                self.purge_directory(methoddir)
                for gene, align in self.alignment[method].items():
                    AlignIO.write(align, os.path.join(methoddir, gene+".fasta") , format=format)

        if self.edited :
            for method in self.edited.keys():
                methoddir = os.path.join(self.aligndir,method)
                self.purge_directory(methoddir)
                for gene, edited in self.edited[method].items():
                    AlignIO.write(edited, os.path.join(methoddir, gene+".fasta") , format=format)


    def concat_align(self):
        """concat aligned sequence"""
        alignment = self.alignment
        filename = "alignstat.html"
        if(self.edited):
            alignment = self.edited
            filename = "editstat.html"
        best_align = self.selection_alignment(alignment, filename=filename)
        if not(best_align):
            print("Could not select best alignment\n")
            return

        if(self.order == []):
            geneorder = sorted(self.genelist)
        else:
            geneorder = self.order

        self.concatfile = raw_input("Enter the final concat alignment filename : ")

        self.concat, self.concatpart = self.concatenateSequences(geneorder, self.speclist, alignment[best_align])


    def save_concat(self, format="phylip-relaxed"):
        """save the concatenated file """
        if not self.concatfile:
            print("Concatfile is not a valid path, setting default name to 'alignment.fasta' before saving")
            self.concatfile = "alignment.fasta"

        # access to here assume that concatfile is defined
        if(self.concat):
            AlignIO.write(self.concat, os.path.join(self.aligndir, self.concatfile), format=format)
        if self.concatpart:
            with open(os.path.join(self.aligndir, self.concatfile.split('.')[0]+'_partition.txt'), 'w') as PART:
                for (genename, pos1, pos2) in self.concatpart:
                    PART.write("LG, %s = %d-%d\n"%(genename, pos1, pos2))


    @classmethod
    def purge_directory(cls, dirname):
        """Empty a repertory"""
        shutil.rmtree(dirname, ignore_errors=True)
        os.makedirs(dirname)


    def edit_data(self):
        """
        Edit sequence, in order to remove genes or species
        """
        print( """

            ---------------------------
            Edit mode
            ---------------------------

            """)

        def validate_choice(inp):
            return inp.lower().startswith('y')

        nooverwrite = True
        if(self.state > 1 or len(self.genelist)<len(self.gene2spec.keys())):
            print("\nIt's seems that you already run edit data on gene once!\nYou can\
             discard your previous data and start over (default) or edit it.")
            nooverwrite = validate_choice(raw_input('Edit previous data y/n? '))

        big_gene_list = self.gene2spec.keys()
        if(nooverwrite):
            big_gene_list = self.genelist

        # part 1 : remove strange gene
        inp = raw_input("1 : Do you want to remove some genes y/n? ")
        to_remove = set([])

        if validate_choice(inp):
            values = raw_input("Percent of gene presence required | list of gene symbol with '-' at the beginning: ")

            while(values and (not ((isFloat(values) and float(values)<=100) or '-' in values))):
                values = raw_input("Not a valid choice, re-enter value or enter to skip : ")
                print values
                values.strip()

            if(isFloat(values)):
                values = float(values)
                required = abs(len(self.spec2genes.keys())*values/100.)
                for gene in big_gene_list :
                    if len(self.gene2spec[gene]) < required:
                        to_remove.add(gene)

            elif ('-' in values):
                gene_regs = [x.strip() for x in values.split('-') if x]
                for greg in gene_regs:
                    greg = re.compile(greg)
                    for g in big_gene_list:
                        if greg.match(g):
                            to_remove.add(g)

            prefix = '-' if to_remove else ''
            print ("The following genes will be removed: " + prefix + "-".join(to_remove))
            answer = raw_input("Please confirm y/n : ")
            if validate_choice(answer):
                self.genelist = list(set(big_gene_list) - to_remove)

        # part 2 remove some species
        to_remove = set()
        inp = raw_input("2 : Do you want to remove some species y/n?")
        if validate_choice(inp):
            values = raw_input("Enter species number separated by '-' : ")
            values = [int(x.strip()) for x in values.split('-') if isInt(x.strip())]
            while(not isinstance(values,list) or len(values)<1):
                values = raw_input("Not a valid choice, re-enter value or enter to skip : ")
                try:
                    values = [int(x.strip()) for x in values.split('-') if isInt(x.strip())]
                except:
                    values = ""

            for ind in values:
                to_remove.add(self.speclist[ind])

            prefix = '-' if to_remove else ''
            print ("The following species will be removed: "+prefix+"-".join(to_remove))
            answer = raw_input("Please confirm y/n : ")
            if validate_choice(answer):
                self.speclist = list(set(self.speclist) - to_remove)


    def sequenceDisplay(self):
        """pretty printing of the gene hits in each genome"""
        #Change names to use only the first elements
        #Get longest species name
        maxInputName = -1
        for s in self.speclist:
            if s:
                if len(s) > maxInputName:
                    maxInputName = len(s)
        if maxInputName < len("Input Name"): maxInputName = len("Input Name")
        maxInputName += 1
        #Setup geneNames lengths
        glengths = [len(x) for x in self.genelist]
        for i,each in enumerate(glengths):
            if each < 6:
                glengths[i] = 6

        print("\nPrinting sequence info.")
        print("\nSequence summary:\n")
        header = "Sp. ID " + "Input name".ljust(maxInputName)
        for i, each in enumerate(self.genelist):
            header += each.ljust(glengths[i])
        print(header)
        for i, key in enumerate(self.speclist):
            row = str(i).ljust(len("Seq ID ")) + str(key).ljust(maxInputName)
            for k, gene in enumerate(self.genelist):
                has_gene = key in self.gene2spec[gene].keys()
                row += str(len(self.gene2spec[gene][key]) if has_gene else 0).ljust(glengths[k])
            print(row)
        print("")


    def report_sequence(self, output=None, show=False):
        """use report lab to save gene_specie map"""

        if not output:
            output = os.path.join(self.outpdf,"matrice.pdf")
        ngenes = len(self.genelist)
        nspecies = len(self.speclist)
        doc = SimpleDocTemplate(output, pagesize=(ngenes*1.65*cm+15*cm, nspecies*cm + 5*cm), rightMargin=5,leftMargin=5, topMargin=10,bottomMargin=5)
        elements = []
        data = []
        data.append([""]+self.genelist)
        for i, key in enumerate(self.speclist):
            row = ["%d - %s"%(i, key)]
            for k, gene in enumerate(self.genelist):
                has_gene = key in self.gene2spec[gene].keys()
                row.append(str(len(self.gene2spec[gene][key]) if has_gene else 0))
            data.append(row)

        #TODO: Get this line right instead of just copying it from the docs
        style = TableStyle([('ALIGN',(1,0),(-1,-1),'RIGHT'),
                               ('TEXTCOLOR',(0,0),(0,-1),colors.red),
                               ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
                               ('TEXTCOLOR',(0,0),(-1,0),colors.blue),
                               ('ALIGN',(0,0),(-0,-1),'CENTER'),
                               ('INNERGRID', (0,0), (-1,-1), 0.1, colors.gray),
                               ('BOX', (0,0), (-1,-1), 0.25, colors.white),
                               ('LINEAFTER', (0,0), (0,-1), 2, colors.black),
                               ('LINEBELOW', (0,0), (-1,0), 2, colors.black),
                               ])

        #Configure style and word wrap
        s = getSampleStyleSheet()
        s = s["BodyText"]
        s.wordWrap = 'CJK'
        data2 = [[Paragraph(cell, s) for cell in row] for row in data]
        t=Table(data2, colWidths=[8*cm]+ngenes*[1.6*cm])
        t.setStyle(style)
        #Send the data and build the file
        elements.append(t)
        doc.build(elements)
        if(show):
            err = executeCMD(pdfviewer + " "+output +" &", "evince", useos=True)
            if (err):
                print(pdfviewer+" not working, Failling back silently"+"\n"+err)


    def perform_align(self):
        print("""
            ---------------------------
            Multiple Alignment mode
            ---------------------------
        """)
        error, binaries = check_binaries([muscle, mafft, clustalo, prank, fsa, probcons], getAll=False)
        if(error):
            print(error)
            print("Alignment not performed!")
            return
        else :
            global ALIGN_METHODS
            ALIGN_METHODS = dict((k,list(set(v).difference(binaries))) for (k,v) in ALIGN_METHODS.items() if k not in binaries)

        print("\nSelect alignment ( " + " - ".join(ALIGN_METHODS.keys()) +" ) or enter to skip \n")

        if not self.alignmethod:
            self.alignmethod = 'muscle'

        error = True
        while error:
            align_prompt = raw_input("Alignment method (default="+self.alignmethod+ "): ")
            align_prompt = align_prompt.strip()
            if align_prompt:
                if align_prompt in ALIGN_METHODS.keys():
                    self.alignmethod = align_prompt
                    error = False
                else:
                    print("Invalid choice "+align_prompt + "! please try again.")
            else:
                print("Default alignment chosen\n")
                error = False

        # if there was any fasta file in the directory
        # it will be deleted
        self.save_genes()

        # return a dict of dict with alignment method as keys and the gene id as
        # the second key
        for method in ALIGN_METHODS[self.alignmethod] :
            if (method != 'hmmmuscle'):
                self.alignment[method] = self.alignSequences(method=method, timeout=99999999, \
                                                                    outdir=None)
            else:
                self.alignment[method] = self.muscleAndHmm(timeout=99999999, outdir=None, loop=hmmiters)
                self.edited[method] = self.alignment[method]
        # state = alignment
        self.state = 2


    def alignSequences(self, method='muscle', outdir=None, timeout=99999999, verbose=True):
        """
        Align the sequence of each gene
        """

        if not outdir:
            outdir =  os.path.join(self.aligndir, method)
            self.purge_directory(outdir)

        output = {}
        alignedSomething = []
        outform = {'muscle': 'fasta' , 'mafft' : 'fasta', 'fsa':'fasta', 'prank':'fasta',\
                'clustalo' : 'fasta', 'probcons' : 'fasta'}
        print("\nStarted Alignment")

        fsa_nice_val = "nice -5"
        prank_nice_val = fsa_nice_val
        probcons_nice_val = ""
        # could have made a loop (lazyness overflow!!)
        # bad manu!
        if 'fsa' in preferences.keys() and preferences['fsa']:
            for s in preferences['fsa']:
                if(s.startswith('nice')):
                    fsa_nice_val = s

        if 'probcons' in preferences.keys() and preferences['probcons']:
            for s in preferences['probcons']:
                if(s.startswith('nice')):
                    probcons_nice_val = s

        if 'prank' in preferences.keys() and preferences['prank']:
            for s in preferences['prank']:
                if(s.startswith('nice')):
                    prank_nice_val = s

        for i, gene in enumerate(self.genelist):
            if verbose: print("... trying to align gene : "+gene+ " ....with "+method)
            geneOutput = None
            inputFile = os.path.join(self.fastadir, gene+ '.fasta')
            outputFile = os.path.join(outdir, gene+"."+outform[method])
            commandLine = ""

            if method == 'muscle':
                commandLine = muscle+' -in ' + inputFile + " -out " + outputFile +" "
                if('muscle' in preferences.keys() and preferences['muscle']):
                    commandLine += " ".join(preferences['muscle'])

            elif method == 'mafft':
                if('mafft' in preferences.keys() and preferences['mafft']):
                    commandLine = mafft + ' ' + " ".join(preferences['mafft'])+ " "
                else:
                    commandLine = mafft+ ' --auto '
                commandLine += inputFile + " > " + outputFile


            if method == 'clustalo':
                commandLine = clustalo + ' -i ' + inputFile + " -o " + outputFile + " -v "
                if('clustalo' in preferences.keys() and preferences['clustalo']):
                    commandLine += " ".join(preferences['clustalo'])

            if method == 'fsa':
                commandLine = fsa_nice_val + " " + fsa + " " + inputFile + " --gui "
                if ('fsa' in preferences.keys() and preferences['fsa']):
                    commandLine += " ".join([x for x in preferences['fsa'] if 'nice' not in x])
                commandLine += " 1> " + outputFile + " 2>/dev/null"
                print commandLine

            if method == 'probcons':
                commandLine = probcons_nice_val + " "  + probcons + " " + inputFile + " "
                if ('probcons' in preferences.keys() and preferences['probcons']):
                    commandLine += " ".join([x for x in preferences['probcons'] if 'nice' not in x])
                commandLine += " > " + outputFile

            if method == 'prank':
                commandLine = prank_nice_val + " " + prank +' -d=' + inputFile + " -o=" + outputFile + " "
                if('prank' in preferences.keys() and preferences['prank']):
                    commandLine += " ".join(preferences['prank'])

            al, res= call_seqprog(commandLine, outputFile, timeout, outputFormat=outform[method])
            output[gene] = res
            alignedSomething.append(al)

        if not any(alignedSomething):
            raise RuntimeError("Nothing was aligned !!!")
        return output

    def get_align_path(self, method, gene):
        """Get the aligned file of a gene"""
        return os.path.join(self.aligndir, method, gene+".fasta")

    def selection_alignment(self, alignment, filename="alignstat.html"):

        """Selection of the best alignment"""
        # what should this method do :
        # verify that the sequence are aligned
        # show alignment statistics
        # ask for user alignment choice
        # remimber user that he should manually check the alignment before
        if self.state < 2 :
            print("Alignment state not reached yet !!")
            return
        print("""

            --------------------------------
            Multiple Alignment Selection
            --------------------------------

            """)
        stats = collections.defaultdict(dict)
        dispo_align = alignment.keys()
        for meth in dispo_align:
            stats[meth]= alignstat(alignment[meth], self.genelist)

        out = PrettyTable(['Genes', 'Align_len', 'Gap_min', 'Gap_max', 'Gap_mean', 'Gap_std', 'Avg ID'])
        out.hrules = ALL
        err, bina = check_binaries([eslalistat])
        for gene in self.genelist:
            try:
                length = []
                gapmin = []
                gapmax = []
                gapmean = []
                gapstd = []
                avgid = []
                for meth in dispo_align:
                    pipe = TerminationPipe(eslalistat+ " " +self.get_align_path(meth, gene), 999, True)
                    pipe.run()
                    if(not (err or pipe.failure)):
                        match = re.search('(?<=Average identity:)\s*[0-9]*\%', pipe.stdout)
                        avgid.append('%s : %s'%(meth, (match.group(0).strip() if match else 'N/A')))
                    else :
                        avgid.append('%s : %s'%(meth, 'N/A'))

                    length.append('%s : %d'%(meth, stats[meth][gene]['length']))
                    gapmin.append('%s : %d'%(meth, stats[meth][gene]['gapmin']))
                    gapmax.append('%s : %d'%(meth, stats[meth][gene]['gapmax']))
                    gapstd.append('%s : %.2f'%(meth, stats[meth][gene]['gapstd']))
                    gapmean.append('%s : %.2f'%(meth, stats[meth][gene]['gapmean']))

                data = [gene, "\n".join(length), "\n".join(gapmin), "\n".join(gapmax),\
                            "\n".join(gapmean), "\n".join(gapstd),"\n".join(avgid)]
                out.add_row(data)
            except KeyError:
                # the gene was not aligned, just pass
                pass

        html = out.get_html_string(attributes = {"class": "pretty"})
        with open(os.path.join(self.outpdf, filename), 'w+') as IN:
            IN.write(html)
        print(out)
        time.sleep(1)

        al_chosen = 'muscle' if 'muscle' in dispo_align else dispo_align[0]
        print('Alignment selection (default = %s) :'%al_chosen)
        print(" -- ".join(dispo_align))
        print("You can run metal between each pair of alignment for comparision\n")
        sel = False
        while(not sel):
            selected = raw_input("Please select an alignment to continue : ")
            selected = selected.strip().lower()
            if(not selected or (selected and selected in dispo_align)):
                sel = True

        return al_chosen if not selected else selected


    def refine_align(self):
        best_align = self.selection_alignment(self.alignment)
        if not(best_align) or 'hmm' in best_align:
            print("Could not select best alignment\nIf you performed the alignment with hmmmuscle, there is nothing to do")
            return
        print("""

            --------------------------------
            Multiple Alignment Edition mode
            --------------------------------
        """)
        methods = {'gblocks':['gblocks'], 'trimal' : ['trimal'], 'both': ['gblocks', 'trimal']}
        error, binaries = check_binaries([gblocks, trimal], getAll=False)
        if(error):
            print(error)
            print("Alignment not refined!")
            return
        else :
            method = dict((k,list(set(v).difference(binaries))) for (k,v) in methods.items() if k not in binaries)

        print("\nSelect editors (" + " - ".join(methods.keys())+ " or enter to skip\n")
        if not self.editmethod:
            self.editmethod = 'gblocks'
        error = True
        while error:
            prompt = raw_input("Alignment method (default="+self.editmethod+ "): ")
            prompt = prompt.strip()
            if prompt:
                if prompt in methods.keys():
                    self.editmethod = prompt
                    error = False
                else:
                    print("Invalid choice "+prompt + "- please try again.")
            else:
                print("Default Editors chosen\n")
                error = False

        for method in methods[self.editmethod]:
            self.edited[method] = self.cleanAlignment(method, best_align, timeout=99999999,\
                                                            outdir=None)
        self.state = 3


    def cleanAlignment(self, method, align, outdir=None, verbose=True, timeout=None):

        if not outdir:
            outdir =  os.path.join(self.aligndir, method)
            self.purge_directory(outdir)
        output = {}
        editedSomething = []
        outform = {'gblocks': 'fasta', 'trimal':'fasta'}
        gblocks_out= []

        for i, gene in enumerate(self.genelist):
            if verbose: print("... trying to align gene : "+gene+ " ....with "+method +"  \r")
            geneOutput = None
            inputFile = os.path.join(self.aligndir, align, gene+ '.fasta')
            outputFile = os.path.join(outdir, gene+"."+outform[method])
            commandLine = ""

            if 'trimal' == method:
                commandLine = trimal + " -in " + inputFile + " -out " + outputFile + " -fasta "
                if('trimal' in preferences.keys() and preferences['trimal']):
                    commandLine += " ".join(preferences['trimal'])
                else:
                    commandLine += '-automated1'


            elif 'gblocks' == method:
                outputFile = inputFile + '-gb'
                commandLine = gblocks + " " + inputFile + " "
                gblocks_out.append(outputFile)
                if('gblocks' in preferences.keys() and preferences['gblocks']):
                    commandLine += " ".join(preferences['gblocks'])

            al, res= call_seqprog(commandLine, outputFile, timeout, outputFormat=outform[method])
            output[gene] = res
            editedSomething.append(al)

        if method == 'gblocks':
            for out in gblocks_out:
                out_basename = os.path.basename(out).replace('-gb', '')
                shutil.move(out, os.path.join(outdir, out_basename))

        if not any(editedSomething):
            raise RuntimeError("Nothing was aligned !!!")
        return output


    @staticmethod
    def concatenateSequences(geneorder, speclist, alignment, missing='-', alpha=generic_protein):
        """Concatenate alignment"""

        align_dict = {}
        genepos = []
        lastpos = 0
        for gene in geneorder:
            if gene in alignment.keys():
                cur_al = alignment[gene]
                al_len = cur_al.get_alignment_length()
                genepos.append((gene, lastpos+1, al_len+lastpos))
                lastpos += al_len
                spec_dict = dict((x.id.split('#')[1].replace('_', ' '), x) for x in cur_al)
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
            align_dict[spec].id = spec.replace(' ', '_')

        return MultipleSeqAlignment(align_dict.values()), genepos


    def muscleAndHmm(self, timeout, outdir, loop=10, minqual=9, check_muscle=True):
        """Align and refine at the same time with hmmalign and muscle"""
        error, binaries = check_binaries([muscle, hmmbuild, hmmalign, eslalimask, eslalimanip])
        if(error):
            print error
            print("Alignment not done!")
            return

        def accessQuality(alignfile, output, minqual):
            # this serve to format the output
            cmd = eslalimanip + " -o %s %s"%(output, alignfile)
            executeCMD(cmd, "alimanip")
            # read output file and extract position
            ppcons = ''
            with open(output, 'r') as text:
                for line in text:
                    if '#=GC' in line and 'PP_cons' in line:
                        ppcons += line.split()[-1].strip()
            highpos = [1 for x in ppcons if x=='*' or (isInt(x) and int(x)>=minqual)]
            return sum(highpos)


        if not outdir:
            outdir =  os.path.join(self.aligndir, 'hmmmuscle')
            self.purge_directory(outdir)

        muscledir = os.path.join(self.aligndir, 'muscle')
        # if we have already done some alignment with muscle
        # we should not waste time
        if (check_muscle):
            if (os.path.exists(muscledir) and os.listdir(muscledir)):
                ans = raw_input('Muscle alignment found\nDo you want to use it (Y/n)? ')
                if ans.strip().upper() != 'Y':
                    self.alignSequences(method='muscle')

        cur_dir = []
        for i in xrange(loop):
            d = os.path.join(outdir, '%d'%i)
            self.purge_directory(d)
            cur_dir.append(d)

        bestiter = {}
        alignment = {}
        for gene in self.genelist:
            inputFile = os.path.join(muscledir, gene+".fasta")
            print('... trying to run hmmbuild and hmmalign '+str(loop)+" times for "+ gene)
            quality = []
            outlist = []
            for i in xrange(loop):
                hmmfile = os.path.join(cur_dir[i],gene+".hmm" )
                outputFile = os.path.join(cur_dir[i],gene )
                buildline = hmmbuild + " %s %s"%(hmmfile, inputFile)
                executeCMD(buildline, 'hmmbuild')
                # will continue if not exception is found
                alignline = hmmalign + " -o %s %s %s"%(outputFile,hmmfile,inputFile)
                executeCMD(alignline, 'hmmalign')
                # finding quality
                eslout = os.path.join(cur_dir[i], gene+".esl")
                quality.append(accessQuality(outputFile, eslout, minqual))
                outlist.append(eslout)
                # next input for hmm is current output
                inputFile = outputFile

            bestiter[gene] = outlist[np.asarray(quality).argmin()]
            cmdline = eslalimask + " "
            align_file = os.path.join(outdir, gene+".sto")
            option = '-o %s -p '%(align_file)
            if 'hmmmuscle' in preferences.keys() and preferences['hmmmuscle']:
                ppcon = False
                for val in preferences['hmmmuscle']:
                    if('ppcons' in val):
                        option += val
                        ppcon = True
                if not ppcon:
                    option += " ".join(val for val in preferences['hmmmuscle'])
            cmdline += option +" "+ bestiter[gene]
            executeCMD(cmdline, "alimask")
            alignment[gene] = AlignIO.read(align_file, format="stockholm")
        return alignment


    def run_raxml(self):
        pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="pyPhylo - An attempt to automatize sequence alignment.", epilog="Help at https://github.com/maclandrol/pyphylo")
    parser.add_argument("--species", "-s", type=argparse.FileType('r'), dest='species', help="Binomial names of species, each on a new line")
    parser.add_argument("-wd", default='Data', dest='workdir', help="Working directory for all output files")
    #parser.add_argument("-alignment", "-a",, dest='align',  help="Alignment method")
    #parser.add_argument("-options", "-o", help="Options file giving detailed instructions to phyloGen.")
    parser.add_argument("-genes", '--fastafile', "-g", dest="pepper", help="The pepper file containing the genes to search for (multiple genes are comma-separated)")
    parser.add_argument("--genome", '-G', dest="genome", help="The masterfile outputed by mfannot.")
    parser.add_argument("--genbank", '--database', '-db', dest="genbank", help="genbank file containing all the data")
    parser.add_argument("--corr", '--corrdb', dest="corrdb", help="Mapping between genome and accession number")
    parser.add_argument("--clades", action="store_true", help="Indicates the 'species' listed in the '-species' file are actually clade name. pyPhylo will try to retrieve every specie in the clade")
    parser.add_argument("--taxonid", action="store_true", dest="taxon", help="Indicates if the specie list is a taxa id list")
    parser.add_argument("--blast", action="store_true", help="blast each result.")
    parser.add_argument("--checkspecie", dest="checkspec", action="store_true", help="blast each result.")
    parser.add_argument("--orthomcl", action="store_true", help="Perform orthomcl clustering to refine data.")
    parser.add_argument("-delay", help="Delay (seconds) when pausing between any internet API calls.")
    parser.add_argument("--concat", help="concatenation order - default : alphabetic order",  type=argparse.FileType('r'))

    args = parser.parse_args()
    specielist = [line.strip() for line in args.species if (line.strip() and not line.startswith('#'))]

    if (args.checkspec):
        specielist = check_specie_list(specielist, args.taxon, args.clades)
    orispec = specielist.pop(0)
    #random.shuffle(specielist)
    genes, genome = parse_master_file(args.genome)
    sequences = SeqIO.parse(args.pepper, format="fasta", alphabet=generic_protein)
    oriseq = {}
    glist = []
    for seq in sequences:
        seq_id = re.search('(?<=ms\s).*(?=\s\;)', seq.description).group(0)
        seq.name = seq_id
        seq.id = seq_id
        seq.seq = seq.seq.strip('*')
        oriseq[seq_id] = seq

    annote_genes(genome, genes, oriseq.keys())
    if args.blast:
        for record in oriseq.values():
            blast_sequence(record, args.workdir)

    concat = []
    if(args.concat):
        concat = [line.strip() for line in args.concat if not line.startswith('#')]
    pyph = pyPhylo(args.workdir, specielist, args.genbank, args.corrdb, (oriseq, orispec), orthomcl=args.orthomcl, concat_order=concat)
