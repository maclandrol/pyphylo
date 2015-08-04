#TODO: use metal to compare alignment
# metal should be disponible for alignment refinment and also after concatenation
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
# If you let this empty, the program will ask you
# about the concatenated alignment filename
concatfile = ""

#binaries
pdfviewer = 'evince'
trimal = 'trimal'
gblocks = 'Gblocks'
raxml = 'raxml'
fsa = 'fsa'
muscle = 'muscle'
prank = 'prank'
clustalo = 'clustalo'
mafft = 'mafft'
probcons = 'probcons'
orthomclpipe = 'orthomcl-pipeline'
hmmbuild = 'hmmbuild'
hmmalign = 'hmmalign'
eslalimask = 'esl-alimask'
eslalimanip = 'esl-alimanip'
# not necessary
eslalistat = 'esl-alistat'
metal ='metal'

hmmiters = 50
#programmes preferences
# aditionnal parameters for each program
# IMPORTANT : do not touch the output format !!
# for fsa, you could add a nice parameter : default will be -5
# for gblocks you should specify the type (-t=p)
# please use lower case for key
# see esl-alimask for hmmmuscle parameter ('--ppcons, --pthresh and --pfract' )
preferences = {
    'gblocks' : ["-t=p", "-b3=5", "-b4=5",  "-b5=h"],
    'trimal' : ['-automated1'],
    'hmmmuscle' : ["--pthresh 0.90", "--pfract 0.5"]
}


#othomcl-pipeline config file
def_conf_file = "/home/manu/Documents/Projects/Programmation/Scripting/pyphylo/orthomcl-pipeline/orthomcl.conf"
