#misc
email = "fmr.noutahi@umontreal.ca"

#directories
dbdir = "dbfiles"
aligndir = "alignments"
fastadir = "fastafile"
ortho_indir = "proteome"
ortho_outdir = "orthomcl"
outputpdf = "outpdf"
outputtree = "outtree"
datadumping = "classdump"

#file
seqdb = "seqs.db"

#binaries
pdfviewer = 'evince'
trimal = 'trimal'
#docs http://trimal.cgenomics.org/use_of_the_command_line_trimal_v1.2
gblocks = 'Gblocks'
#docs
raxml = 'raxml'
metal ='metal'
hmmbuild = 'hmmbuild'
hmmalign = 'hmmalign'
fas = 'fas'
muscle = 'muscle'
prank = 'prank'
clustalo = 'clustalo'
mafft = 'mafft'
orthomclpipe = 'orthomcl-pipeline'

#programmes preferences
# aditionnal parameters for each program
# do not touch the output format !!
# for fsa, you could add a nice parameter : default will be -5
# for gblocks you should specify the type (-t=p)
# please use lower case for key
preferences = {
    'gblocks' : ["-t=p", "-b3=5", "-b4=5",  "-b5=h"],
    'trimal' : ['-automated1']
}

#othomcl-pipeline config file
def_conf_file = "/home/manu/Documents/Projects/Programmation/Scripting/pyphylo/orthomcl-pipeline/orthomcl.conf"
